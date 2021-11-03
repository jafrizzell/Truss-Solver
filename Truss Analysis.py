import numpy as np
import copy


class Structure:
    def __init__(self, dof):
        self.matrix = np.zeros((dof, dof))
        self.forces = []
        self.displacements = []
        self.members = []

    def add_member(self, member):
        self.members.append(member)
        if len(member.dof) == 4:
            for i in range(len(member.k_matrix)):
                for j in range(len(member.k_matrix[i])):
                    row = member.dof[i]
                    col = member.dof[j]
                    self.matrix[row-1][col-1] += round(member.ea*member.k_matrix[i][j]/member.l, 5)
            for row in range(len(member.k_matrix)):
                print(member.ea*np.asarray(member.k_matrix[row])/member.l, member.dof[row])
            print('')
        if len(member.dof) == 6:
            for i in range(len(member.k_matrix)):
                for j in range(len(member.k_matrix[i])):
                    row = member.dof[i]
                    col = member.dof[j]
                    self.matrix[row-1][col-1] += round(member.e*member.k_matrix[i][j]/member.l, 5)
            for row in range(len(member.k_matrix)):
                print(member.e*np.asarray(member.k_matrix[row])/member.l, member.dof[row])
            print('')

    def solve(self):
        disp_copy = copy.deepcopy(self.displacements)
        force_copy = copy.deepcopy(self.forces)
        ind_unknown = []
        dof_unknown = []
        count = 0
        for a in range(len(self.displacements)):
            if type(self.displacements[a]) == str:
                ind_unknown.append(a)
                dof_unknown.append(a+1)

        num_unknown = len(ind_unknown)
        num_known = len(self.displacements) - num_unknown
        for i in ind_unknown:
            self.matrix[[i, count]] = self.matrix[[count, i]]
            self.matrix[:, [i, count]] = self.matrix[:, [count, i]]
            self.displacements[i], self.displacements[count] = self.displacements[count], self.displacements[i]
            self.forces[i], self.forces[count] = self.forces[count], self.forces[i]
            count += 1

        print(np.asarray(self.matrix))
        k11 = self.matrix[:num_unknown, :num_unknown]
        k12 = self.matrix[:num_unknown, -num_known:]
        k21 = self.matrix[num_unknown:, :num_unknown]
        k22 = self.matrix[-num_known:, -num_known:]


        print("\n K11")
        print(np.asarray(k11))
        print("\n K12")
        print(np.asarray(k12))
        print("\n K21")
        print(np.asarray(k21))
        print("\n K22")
        print(np.asarray(k22))
        print('')

        dk = self.displacements[num_unknown:]
        qk = self.forces[:num_unknown]
        du = np.dot(np.linalg.inv(k11), (qk-np.dot(k12, dk)))

        qu = np.dot(k21, du) + np.dot(k22, dk)
        qu = qu.tolist()

        for f in range(len(ind_unknown)):
            for i in range(len(ind_unknown)):
                disp_copy[ind_unknown[i]] = du[i]
            qu.insert(ind_unknown[f], force_copy[ind_unknown[f]])

        for num in range(len(qu)):
            disp_copy[num], qu[num] = round(disp_copy[num], 5), round(qu[num], 5)

        print('Displacement =', du, 'on directions', dof_unknown, '\n')
        print("System Displacements:", disp_copy, '\n')
        print("Nodal Reactions:", qu, '\n')

        mem_num = 1
        print("Member internal forces (positive is compression, negative is tension)\n")
        for member in self.members:
            # TODO: Add proper support for the forces inside a beam element with f = k_element * d_element
            if len(member.k_matrix) == 4:
                du = [disp_copy[member.nx-1], disp_copy[member.ny-1], disp_copy[member.fx-1], disp_copy[member.fy-1]]
                member_mat = [member.lam_x, member.lam_y, -member.lam_x, -member.lam_y]
                member_force = round(member.ea/member.l * np.dot(member_mat, du), 5)
            if len(member.k_matrix) == 6:
                member_force = 'n/a'
            else:
                du = [disp_copy[member.nx-1], disp_copy[member.fx-1]]
                member_force = member.ea / member.l * (du[1] - du[0])
            print('Force in member #', mem_num, 'is', member_force, '\n')
            mem_num += 1


class Member:
    class Spring:
        def __init__(self):
            self.nx = 0
            self.fx = 0
            self.ea = 0
            self.l = 0
            self.lam_x = 0
            self.lam_y = 0
            self.k_matrix = np.zeros((2,2))
            self.dof = [0, 0]

        def define(self, nx, fx, k, l, lam_x, lam_y):
            self.nx = nx
            self.fx = fx
            self.ea = k
            self.l = l
            self.lam_x = lam_x
            self.lam_y = lam_y
            self.dof = [nx, fx]

            self.k_matrix = [
                [1, -1],
                [-1, 1]
            ]

    class Truss:
        def __init__(self):
            self.nx = 0
            self.ny = 0
            self.fx = 0
            self.fy = 0
            self.ea = 0
            self.l = 0
            self.lam_x = 0
            self.lam_y = 0
            self.k_matrix = np.zeros((4, 4))
            self.dof = [0, 0, 0, 0]

        def define(self, nx, ny, fx, fy, ea, l, lam_x, lam_y):
            self.nx = nx
            self.ny = ny
            self.fx = fx
            self.fy = fy
            self.ea = ea
            self.l = l
            self.lam_x = lam_x
            self.lam_y = lam_y
            self.dof = [nx, ny, fx, fy]

            self.k_matrix = \
            [
                [self.lam_x**2, self.lam_x*self.lam_y, -1*self.lam_x**2, -1*self.lam_x*self.lam_y],
                [self.lam_x*self.lam_y, self.lam_y**2, -1*self.lam_x*self.lam_y, -1*self.lam_y**2],
                [-1*self.lam_x**2, -1*self.lam_x*self.lam_y, self.lam_x**2, self.lam_x*self.lam_y],
                [-1*self.lam_x*self.lam_y, -1*self.lam_y**2, self.lam_x*self.lam_y, self.lam_y**2]
            ]

    class Beam:
        def __init__(self):
            self.nx = 0
            self.ny = 0
            self.nz = 0
            self.fx = 0
            self.fy = 0
            self.fz = 0
            self.e = 0
            self.a = 0
            self.i = 0
            self.l = 0
            self.lam_x = 0
            self.lam_y = 0
            self.k_matrix = np.zeros((6, 6))
            self.dof = [0, 0, 0, 0, 0, 0]

        def define(self, nx, ny, nz, fx, fy, fz, e, a, i, l, lam_x, lam_y):
            self.nx = nx
            self.ny = ny
            self.nz = nz
            self.fx = fx
            self.fy = fy
            self.fz = fz
            self.e = e
            self.a = a
            self.i = i
            self.l = l
            self.lam_x = lam_x
            self.lam_y = lam_y
            self.dof = [nx, ny, nz, fx, fy, fz]
            self.k_matrix = \
                [
                    [self.a*self.lam_x**2 + 12*self.i/self.l**2*self.lam_y**2, (self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, -6*self.i/self.l*self.lam_y, -(self.a*self.lam_x**2+12*self.i/self.l**2*self.lam_y**2), -(self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, -6*self.i/self.l*self.lam_y],
                    [(self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, self.a*self.lam_y**2+12*self.i/self.l**2*self.lam_x**2, 6*self.i/self.l*self.lam_x, -(self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, -(self.a*self.lam_y**2+12*self.i/self.l**2*self.lam_x**2), 6*self.i/self.l*self.lam_x],
                    [-6*self.i/self.l*self.lam_y, 6*self.i/self.l*self.lam_x, 4*self.i, 6*self.i/self.l*self.lam_y, -6*self.i/self.l*self.lam_x, 2*self.i],
                    [-(self.a*self.lam_x**2+12*self.i/self.l**2*self.lam_y**2), -(self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, 6*self.i/self.l*self.lam_y, self.a*self.lam_x**2+12*self.i/self.l**2*self.lam_y**2, (self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, 6*self.i/self.l*self.lam_y],
                    [-(self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, -(self.a*self.lam_y**2+12*self.i/self.l**2*self.lam_x**2), -6*self.i/self.l*self.lam_x, (self.a-12*self.i/self.l**2)*self.lam_x*self.lam_y, self.a*self.lam_y**2+12*self.i/self.l**2*self.lam_x**2, -6*self.i/self.l*self.lam_x],
                    [-6*self.i/self.l*self.lam_y, 6*self.i/self.l*self.lam_x, 2*self.i, 6*self.i/self.l*self.lam_y, -6*self.i/self.l*self.lam_x, 4*self.i]
                ]


dof = 12
structure = Structure(dof)

b1 = Member.Beam()
b2 = Member.Beam()
b3 = Member.Beam()

b1.define(1, 2, 3, 4, 5, 6, 210*10**9, 1*10**-2, 2*10**-4, 3, 1, 0)
b2.define(4, 5, 6, 7, 8, 9, 210*10**9, 1*10**-2, 2*10**-4, 4.24, 0.707, -0.707)
b3.define(7, 8, 9, 10, 11, 12, 210*10**9, 1*10**-2, 2*10**-4, 3, 1, 0)

structure.add_member(b1)
structure.add_member(b2)
structure.add_member(b3)

structure.forces = ['f', 'f', 'f', 0, -10000, -5000, 0, -10000, 5000, 'f', 'f', 'f']
structure.displacements = [0, 0, 0, 'f', 'f', 'f', 'f', 'f', 'f', 0, 0, 0]


structure.solve()
