import os, sys
import math
import importlib
import numpy as np
import networkx as nx
import argparse
import mdtraj as md
import json
from multiprocessing import Pool


class protein:
    def __init__(self, reaction, output_dir):
        self.reaction = reaction
        self.output_dir = output_dir
        self.bat = {2: "bond", 3: "angle", 4: "torsion"}

    def load(self):
        self.unique = [tuple([eval(items) for items in line.split()]) for line in open(os.path.join(self.output_dir, "topology.txt")).readlines()]
        self.S1D = np.load(os.path.join(self.output_dir, "S1D.npy"))
        topology = md.load(os.path.join(self.output_dir, "topology.pdb")).top
        self.atom_num = topology.n_atoms
        self.residue_dict = {
            f"{residue}{r + 1}": topology.select(f"resid {r}")
            for r, residue in enumerate(topology.residues)
        }
        self.entropy_map = np.zeros((topology.n_residues, topology.n_residues))
        self.MI_bank = {}

    def get_topology(self, topology):
        graph_dict = {}
        for atom_pair in topology.bonds:
            a1, a2 = atom_pair
            index1 = a1.index
            index2 = a2.index
            if index1 not in graph_dict:
                graph_dict[index1] = []
            if index2 not in graph_dict:
                graph_dict[index2] = []
            graph_dict[index1].append(index2)
            graph_dict[index2].append(index1)

        G = nx.Graph(graph_dict)

        data = [items for items in nx.all_pairs_shortest_path(G, cutoff=3)]

        total_paths = [value for key, dic in data for value in dic.values() if len(value) > 1]

        unique = []
        for items in total_paths:
            if tuple(items) not in unique and tuple(reversed(items)) not in unique:
                unique.append(tuple(items))
        unique.sort()
        unique.sort(key=len)
        with open(os.path.join(self.output_dir, 'topology.txt'), 'w') as fw:
            for items in unique:
                fw.write((" ").join(str(x) for x in items) + "\n")
        return unique

    def dof2bat(self, arg):
        i, X = arg
        save = os.path.join(self.output_dir, "dof%d.npy" % i)
        if len(X) == 2:
            np.save(
                save, (md.compute_distances(self.trajectory, [X]) * 10.0).flatten()
            )  # nm2ang
        elif len(X) == 3:
            np.save(save, md.compute_angles(self.trajectory, [X]).flatten())
        elif len(X) == 4:
            np.save(save, md.compute_dihedrals(self.trajectory, [X]).flatten())

    def xyz2bat(self):
        self.trajectory = md.load("./dataset/%s.pdb" % self.reaction)
        topology = self.trajectory.top
        self.trajectory[0].save_pdb(os.path.join(self.output_dir, "topology.pdb"))
        unique = self.get_topology(topology)
        with Pool() as p:
            p.map(self.dof2bat, [(i, items) for i, items in enumerate(unique)])
            S1D = p.map(
                self.entropy_1D,
                [(i, self.bat[len(unique[i])]) for i in range(len(unique))],
            )
            p.close()
            p.join()
        np.save(os.path.join(self.output_dir, "S1D.npy"), np.array(S1D))

    def get_jacobian(self, X_bin_edge, xjtype, Y_bin_edge, yjtype):
        if Y_bin_edge is None and yjtype is None:
            if xjtype == "bond":
                return X_bin_edge**2
            elif xjtype == "torsion":
                return 1
            else:
                return math.sin(X_bin_edge)
        else:
            return_dict = {
                ("bond", "bond"): X_bin_edge**2 * Y_bin_edge**2,
                ("bond", "angle"): X_bin_edge**2 * math.sin(Y_bin_edge),
                ("bond", "torsion"): X_bin_edge**2,
                ("angle", "angle"): math.sin(X_bin_edge) * math.sin(Y_bin_edge),
                ("angle", "torsion"): math.sin(X_bin_edge),
                ("torsion", "bond"): X_bin_edge**2,
                ("torsion", "angle"): math.sin(X_bin_edge),
            }
            return return_dict.get((xjtype, yjtype), 1)

    def entropy_1D(self, args):
        counts, bin_edges = np.histogram(
            np.load(os.path.join(self.output_dir, "dof%d.npy" % args[0])), 50
        )
        jtype = args[1]
        sample_size = np.sum(counts)
        prob_den = counts / sample_size
        dx = bin_edges[1] - bin_edges[0]
        bin_edges += dx / 2

        entropy_sum = -np.sum(
            [
                prob_den[i]
                * math.log(
                    prob_den[i]
                    / (self.get_jacobian(bin_edges[i], jtype, None, None) * dx)
                )
                for i in range(len(prob_den))
                if prob_den[i] != 0
            ]
        )

        return entropy_sum + (np.count_nonzero(prob_den) - 1) / (2 * sample_size)

    def Mutual_Information(self, args):
        key = args[0]
        if key in self.MI_bank:
            return self.MI_bank[key]
        else:
            H, X_bin_edges, Y_bin_edges = np.histogram2d(
                np.load(os.path.join(self.output_dir, "dof%d.npy" % key[0])),
                np.load(os.path.join(self.output_dir, "dof%d.npy" % key[1])),
                bins=50,
            )
            xjtype, yjtype = args[1], args[2]
            dx, dy = X_bin_edges[1] - X_bin_edges[0], Y_bin_edges[1] - Y_bin_edges[0]
            X_bin_edges += dx / 2
            Y_bin_edges += dy / 2
            sample_size = np.sum(H)
            H = H / sample_size

            entropy_sum = -np.sum(
                [
                    H[row][col]
                    * math.log(
                        H[row][col]
                        / (
                            self.get_jacobian(
                                X_bin_edges[row], xjtype, Y_bin_edges[col], yjtype
                            )
                            * dx
                            * dy
                        )
                    )
                    for row in range(50)
                    for col in range(50)
                    if H[row][col] != 0
                ]
            )

            MI = (
                self.S1D[key[0]]
                + self.S1D[key[1]]
                - (entropy_sum + (np.count_nonzero(H) - 1) / (2 * sample_size))
            )
            self.MI_bank[key] = MI
            return MI

    def get_sorted_indices(self, entropy_list, num_indices):
        return sorted(range(len(entropy_list)), key=lambda i: entropy_list[i])[
            -num_indices:
        ]

    def residue_entropy_prep(self, residue_atoms):
        atom_num = len(residue_atoms)

        bond_indices = [
            index
            for index, tup in enumerate(self.unique)
            if all(elem in residue_atoms for elem in tup) and len(tup) == 2
        ]
        angle_indices = [
            index
            for index, tup in enumerate(self.unique)
            if all(elem in residue_atoms for elem in tup) and len(tup) == 3
        ]
        torsion_indices = [
            index
            for index, tup in enumerate(self.unique)
            if all(elem in residue_atoms for elem in tup) and len(tup) == 4
        ]
        bond_entropy = [self.S1D[i] for i in bond_indices]
        angle_entropy = [self.S1D[i] for i in angle_indices]
        torsion_entropy = [self.S1D[i] for i in torsion_indices]

        bond_mea_indices = self.get_sorted_indices(bond_entropy, atom_num - 1)
        angle_mea_indices = self.get_sorted_indices(angle_entropy, atom_num - 2)
        torsion_mea_indices = self.get_sorted_indices(torsion_entropy, atom_num - 3)
        indices = (
            [bond_indices[i] for i in bond_mea_indices]
            + [angle_indices[i] for i in angle_mea_indices]
            + [torsion_indices[i] for i in torsion_mea_indices]
        )
        return sum([self.S1D[i] for i in indices]), [
            (indices[i], indices[j])
            for i in range(len(indices) - 1)
            for j in range(i + 1, len(indices))
        ], len(indices)

    def residueS(self, resikey):
        residue_atoms = self.residue_dict[resikey]
        return self.residue_entropy_prep(residue_atoms)

    def residue2residue(self, resikey1, resikey2):
        residue_atoms = np.concatenate(
            (self.residue_dict[resikey1], self.residue_dict[resikey2])
        )
        return self.residue_entropy_prep(residue_atoms)

    def mapping(self, temperature):
        conversion_factor = -1.987204259e-3 * temperature  # kcal/mol
        rk = [key for key in self.residue_dict]
        total_S1D, indices, idx_sizes, num_dof = [], [], [], []
        for i in range(len(rk)):
            for j in range(i, len(rk)):
                if i == j:
                    arg1, arg2, arg3 = self.residueS(rk[i])
                elif i != j:
                    arg1, arg2, arg3 = self.residue2residue(rk[i], rk[j])
                total_S1D.append(arg1)
                indices = indices + arg2
                idx_sizes.append(len(arg2))
                num_dof.append(arg3)
        with Pool() as p:
            MI = p.map(
                self.Mutual_Information,
                [
                    (
                        indices[i],
                        self.bat[len(self.unique[indices[i][0]])],
                        self.bat[len(self.unique[indices[i][1]])],
                    )
                    for i in range(len(indices))
                ],
            )
            p.close()
            p.join()
        MIST = []
        pre = 0
        for i, items in enumerate(idx_sizes):
            MI_subset = np.array(MI[pre : pre + items])
            pre = pre + items
            N = num_dof[i]
            matrix = np.zeros((N, N))
            row, col = np.triu_indices(N, k=1)
            matrix[row, col] = MI_subset
            MIST.append(-sum([max(row) for row in matrix]))

        idx = 0
        for i in range(len(rk)):
            for j in range(i, len(rk)):
                if i == j:
                    self.entropy_map[i, j] = conversion_factor * (
                        total_S1D[idx] + MIST[idx]
                    )
                else:
                    self.entropy_map[i, j] = conversion_factor * MIST[idx]
                    self.entropy_map[j, i] = conversion_factor * MIST[idx]
                idx += 1
        np.save(
            os.path.join(self.output_dir, "%s_entropy_map.npy" % self.reaction),
            self.entropy_map,
        )
        for i in range(len(self.S1D)):
            os.remove(os.path.join(self.output_dir, "dof%d.npy" % i))

    def entire(self, temperature):
        conversion_factor = -1.987204259e-3 * temperature  # kcal/mol
        bond_indices = [index for index, tup in enumerate(self.unique) if len(tup) == 2]
        angle_indices = [
            index for index, tup in enumerate(self.unique) if len(tup) == 3
        ]
        torsion_indices = [
            index for index, tup in enumerate(self.unique) if len(tup) == 4
        ]

        bond_entropy = [self.S1D[i] for i in bond_indices]
        angle_entropy = [self.S1D[i] for i in angle_indices]
        torsion_entropy = [self.S1D[i] for i in torsion_indices]

        bond_mea_indices = self.get_sorted_indices(bond_entropy, self.atom_num - 1)
        angle_mea_indices = self.get_sorted_indices(angle_entropy, self.atom_num - 2)
        torsion_mea_indices = self.get_sorted_indices(torsion_entropy, self.atom_num - 3)
        indices = (
            [bond_indices[i] for i in bond_mea_indices]
            + [angle_indices[i] for i in angle_mea_indices]
            + [torsion_indices[i] for i in torsion_mea_indices]
        )

        indices_2D = [(indices[i], indices[j]) for i in range(len(indices) - 1) for j in range(i + 1, len(indices))]

        with Pool() as p:
            MI = p.map(
                self.Mutual_Information,
                [
                    (
                        indices_2D[i],
                        self.bat[len(self.unique[indices_2D[i][0]])],
                        self.bat[len(self.unique[indices_2D[i][1]])],
                    )
                    for i in range(len(indices_2D))
                ],
            )
            p.close()
            p.join()
        N = len(indices)
        matrix = np.zeros((N, N))
        row, col = np.triu_indices(N, k=1)
        matrix[row, col] = MI
        entropy1D_total = sum([self.S1D[i] for i in indices])
        MIE = np.sum(MI)
        MIST = np.sum([max(items) for items in matrix])
        with open(os.path.join(self.output_dir, "entropy.log"), "w") as fw:
            for i in indices:
                fw.write(
                    "%s %s\t %f\n"
                    % (
                        self.bat[len(self.unique[i])],
                        (" ").join(str(x) for x in self.unique[i]),
                        self.S1D[i],
                    )
                )
            for i in range(len(indices) - 1):
                for j in range(i + 1, len(indices)):
                    fw.write(
                        "%s %s\t%s %s\t%f\n"
                        % (
                            self.bat[len(self.unique[indices[i]])],
                            (" ").join(str(x) for x in self.unique[indices[i]]),
                            self.bat[len(self.unique[indices[j]])],
                            (" ").join(str(x) for x in self.unique[indices[j]]),
                            matrix[i, j],
                        )
                    )
            fw.write("MIE : %f\n" % MIE)
            fw.write("MIST : %f\n" % MIST)
            fw.write("MIE Entropy : %f\n" % (entropy1D_total - MIE))
            fw.write("MIST Entropy : %f\n" % (entropy1D_total - MIST))
            fw.write(
                "MIE Entropy (kcal/mol) : %f\n"
                % (conversion_factor * (entropy1D_total - MIE))
            )
            fw.write(
                "MIST Entropy (kcal/mol) : %f\n"
                % (conversion_factor * (entropy1D_total - MIST))
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("--reaction", type=str, default="", help="Name of the reaction")
    parser.add_argument(
        "--temperature", type=float, default=298.15, help="Temperature for EPS"
    )
    parser.add_argument("--job", type=int, default=0, help="job")

    args = parser.parse_args()
    if args.reaction == "":
        print("Error: Reaction Name is not provided.")
        sys.exit()

    output_dir = "./%s" % args.reaction
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    PE = protein(args.reaction, output_dir)
    if args.job == 0:
        PE.xyz2bat()
    elif args.job == 1:
        PE.load()
        PE.entire(args.temperature)
        PE.mapping(args.temperature)
