import numpy as np
import scipy.constants as const
import copy
from src import system, potential1D


class ConveyorBelt:
    '''
    Conveyor belt ensemble class
    organizes the replicas and their coupling
    '''
    def traj_clear(self):
        '''
        deletes trajectories of replicas
        :return: None
        '''
        self.systrajs = []

    def add_sys(self):
        '''
        adds a replica to the ensemble
        :return: None
        '''
        self.systems.append(copy.deepcopy(self.systems[0]))

    def init_mem(self):
        '''
        initializes memory
        :return: None
        '''
        self.num_gp = 11
        self.mem_fc = 0.0001
        #        self.mem=np.array([2.2991 ,  2.00274,  1.84395,  1.83953,  2.0147])
        #        memory for perturbed hosc, alpha=10.0, gamma=0.0, 8 replica, num_gp=6, fc=0.00001, 1E6 steps
        self.mem = np.zeros(self.num_gp - 1)
        self.gp_spacing = self.dis / float(self.num_gp - 1.0)
        self.biasene = 0.0
        # print('Distance: ', self.dis, self.dis / np.pi)
        # print('GP Distance: ', self.gp_spacing, self.gp_spacing / np.pi)
        # print('Gridpoints: ', np.linspace(0, self.num_gp - 1, self.num_gp) * self.gp_spacing)
        # print('Gridpoints: ', np.linspace(0, self.num_gp - 1, self.num_gp) * self.gp_spacing / np.pi)

    def build_mem(self):
        '''
        increments biasing memory
        :return: None
        '''
        active_gp = int(np.floor((self.caplam % self.dis) / self.gp_spacing + 0.5))
        self.mem[active_gp % (self.num_gp - 1)] += self.mem_fc

    def apply_mem(self):
        active_gp = int(np.floor((self.caplam % self.dis) / self.gp_spacing + 0.5))
        dg = (self.caplam % self.dis) / self.gp_spacing - float(active_gp)
        if dg < 0:
            self.biasene = self.mem[(active_gp - 1) % (self.num_gp - 1)] * self.spline(1.0 + dg) + self.mem[
                active_gp % (self.num_gp - 1)] * self.spline(-dg)
        else:
            self.biasene = self.mem[active_gp % (self.num_gp - 1)] * self.spline(dg) + self.mem[
                (active_gp + 1) % (self.num_gp - 1)] * self.spline(1.0 - dg)
        # print("{:5.2f}{:5.2f}{:8.3f}{:3d}{:8.3f}{:8.3f}{:8.3f} {:s}".format(self.caplam, (self.caplam%self.dis),
        # (self.caplam%self.dis)/self.gp_spacing, active_gp,
        # self.gp_spacing*active_gp, dg, ene, np.array2string(self.mem)))

    @staticmethod
    def spline(dg):
        '''
        calculates the value of the spline function depending on the deviation dg from the grid point
        :param dg: deviation from gridpoint (absolute value)
        :return: value of spline (float)
        '''
        if dg < 0.0:
            print('distance smaller than 0')
        elif dg < 1.0:
            return 1.0 - 3.0 * dg * dg + 2 * dg * dg * dg
        else:
            return 0.0

    def calc_lam(self, caplam, i, states=None):
        '''
        calculates lam_i for replica i depending on ensemble state caplam
        :param caplam: state of ensemble (capital lambda) 0 <= caplam < 2 pi
        :param i: index of replica
        :return: lam_i
        '''
        ome = (caplam + i * self.dis) % (2. * np.pi)
        if ome > np.pi:
            ome = 2.0 * np.pi - ome
        return ome / np.pi


    def updateBlam(self, caplam):
        '''
        updates the state of the ensemble and the replicas accordingly
        :param caplam: capital lambda 0 <= caplam < 2 pi
        :return: caplam
        '''
        self.caplam = caplam
        for i in range(self.num):
            self.systems[i].updateLam(self.calc_lam(caplam, i))
        self.apply_mem()
        self.ene = self.calc_ene()
        return caplam

    def calc_ene(self):
        '''
        calculates energy of ensemble
        :return: None
        '''
        ene = 0.0
        for i in range(self.num):
            ene += self.systems[i]._currentTotPot+self.systems[i]._currentTotKin
        self.ene = ene + self.biasene
        return self.ene

    def propagate(self):
        '''
        propagates ensemble
        :return: None
        '''
        self.stepcount += 1
        #        self.oldstate=self.state
        self.state = []
        for j in range(self.num):
            #            self.systems[j].oldpos=self.systems[j].pos
            #            self.systems[j].pos += self.systems[j].randomShift()
            #            self.systems[j].updateEne()
            self.systems[j].propagate()
        oldEne = self.calc_ene()
        oldBiasene = self.biasene
        oldBlam = self.caplam
        self.caplam += (np.random.rand() * 2.0 - 1.0) * np.pi / 4.0
        self.caplam = self.caplam % (2.0 * np.pi)
        self.updateBlam(self.caplam)
        newEne = self.ene
        if newEne < oldEne:
            for i in range(self.num):
                self.state.append(self.systems[i].getCurrentState())
            self.systrajs.append(self.state)
            self.traj.append(np.array([self.stepcount, self.caplam, newEne, self.biasene]))
        elif np.random.rand() <= np.exp(-self.beta * (newEne - oldEne)):
            for i in range(self.num):
                self.state.append(self.systems[i].getCurrentState())
            self.systrajs.append(self.state)
            self.traj.append(np.array([self.stepcount, self.caplam, newEne, self.biasene]))
        else:
            self.reject += 1
            self.updateBlam(oldBlam)
            # self.revert()
            for i in range(self.num):
                self.state.append(self.systems[i].getCurrentState())
            self.systrajs.append(self.state)
            self.traj.append(np.array([self.stepcount, oldBlam, oldEne, oldBiasene]))
        if self.build:
            self.build_mem()

    def revert(self):
        '''
        reverts last propagation step
        :return: None
        '''
        for j in range(self.num):
            self.systems[j].revert()
        self.calc_ene()

    def printTraj(self, fileString, append='ab'):
        '''
        prints the ensemble and replica trajectories to file
        :param fileString: base string for file name
        :param append: append?
        :return: None
        '''
        outfile = open(fileString + '_ens.dat', append)
        np.savetxt(outfile, np.array(self.traj), fmt='%25.15e')
        outfile.close()
        self.traj = []
        self.systrajs = np.array(self.systrajs)
        for i in range(self.num):
            outfile = open(fileString + '_' + str(i) + '.dat', append)
            np.savetxt(outfile, self.systrajs[:, i, :], fmt='%25.15e')
            outfile.close()
        self.systrajs = []
        return

    def initTraj(self, fileString):
        '''
        initializes the output files
        :param fileString: base string for output file names
        :return: None
        '''
        outfile = open(fileString + '_ens.dat', 'wb')
        outfile.close()
        for i in range(self.num):
            outfile = open(fileString + '_' + str(i) + '.dat', 'wb')
            outfile.close()
        return

    def print_systems(self):
        '''
        prints state of systems
        :return: None
        '''
        print(self.__str__())
        return

    def __str__(self):
        '''
        :return: ensemble state string
        '''
        outstr = ''
        for i in range(self.num):
            outstr += '{:d}{:10.2f}{:10.3f}\n'.format(i, self.systems[i]._currentLam, self.systems[i].totene)
        return outstr

    def __repr__(self):
        '''
        :return: ensemble state string
        '''
        return self.__str__()

    import src.integrator as integrator

    def __init__(self, caplam, num,
                 system=system.perturbedSystem(temperature=300.0, lam=0.0, potential=potential1D.pertHarmonicOsc1D(fc=1.0, alpha=10.0),
                                               integrator=integrator.metropolisMonteCarloIntegrator()), build=False):
        '''
        initialize Ensemble object
        :param caplam: state of ensemble, 0 <= caplam < pi
        :param num: number of replicas
        :param system: a system1D instance
        :param build: build memory?
        '''
        assert 0.0 <= caplam <= 2 * np.pi, "caplam not allowed"
        assert num >= 1, "At least one system is needed"
        self.num = num
        self.caplam = caplam
        self.stepcount = 0
        self.reject = 0
        self.beta = 1.0 / (const.gas_constant / 1000.0 * system.temperature)
        self.dis = 2.0 * np.pi / num
        self.build = build
        self.systems = [system]
        self.systrajs = []
        self.traj = []
        self.state = [[[]]] * self.num
        #        self.state.append(self.systems[0].getState())
        self.oldstate = []
        for i in range(self.num - 1):
            self.add_sys()

        # initialize memory variables
        self.num_gp = None
        self.mem_fc = None
        self.mem = None
        self.gp_spacing = None
        self.biasene = None
        self.init_mem()

        self.updateBlam(self.caplam)
        while self.stepcount == self.reject:
            self.propagate()
        self.systrajs = []
        self.traj = []
        self.stepcount = -1
        self.reject = 0
        self.ene = 0.0


def calc_traj(steps=1, ens=ConveyorBelt(0.0, 8)):
    '''
    function to propagate the ensemble ens steps steps
    :param steps: (int) steps
    :param ens: Ensemble object
    :return: tuple of replica trajectories and ensemble trajectory
    '''
    for i in range(steps):
        ens.propagate()
    print("Rejected ", ens.reject)
    return np.array(ens.systrajs), np.array(ens.traj)


def calc_traj_file(steps=1, ens=ConveyorBelt(0.0, 8), filestring='traj'):
    '''
    function to propagate the ensemble ens steps steps and write to file with name filestring
    :param filestring: file name string
    :param steps: (int) steps
    :param ens: Ensemble object
    :return: number of rejected steps
    '''
    ens.initTraj(filestring)
    for i in range(steps):
        ens.propagate()
        ens.printTraj(filestring, append='ab')
    return ens.reject
