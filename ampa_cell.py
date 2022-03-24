#!/usr/bin/env python3

import arbor
class single_cell_recipe (arbor.recipe):
    def __init__(self, mech_type, a_ratio=0, b_ratio=0, c_ratio=0):
        arbor.recipe.__init__(self)
        self.props = arbor.neuron_cable_properties()
        self.type = mech_type 

        if self.type == 'approx':
            self.props.catalogue = arbor.load_catalogue("/home/abiakarn/git/synapses/Ampa-catalogue-approx.so")
        elif self.type == 'gt':
            self.props.catalogue = arbor.load_catalogue("/home/abiakarn/git/synapses/Ampa-catalogue-gt.so")
        else: 
            raise ValueError('only approx and gt mech_types allowed')

        self.a_ratio = a_ratio
        self.b_ratio = b_ratio
        self.c_ratio = c_ratio

    def num_cells(self):
        return 1

    def cell_description(self, gid):
        tree = arbor.segment_tree()
        tree.append(arbor.mnpos, arbor.mpoint(0, 0, 0, 2), arbor.mpoint(0, 0, 4, 2), tag=1)

        labels = arbor.label_dict({
            'synapse_site': '(location 0 0.5)', 
            'root': '(root)',
        })

        params = {
            'a_ratio': self.a_ratio,
            'b_ratio': self.b_ratio,
            'c_ratio': self.c_ratio,
        } 

        decor = arbor.decor()

        if self.type == 'approx':
            decor.place('"synapse_site"', arbor.synapse('Ampa', params), 'syn')
        elif self.type == 'gt':
            decor.place('"synapse_site"', arbor.synapse('Ampa'), 'syn')

        decor.discretization(arbor.cv_policy_single())

        return arbor.cable_cell(tree, labels, decor)

    def cell_kind(self, gid):
        return arbor.cell_kind.cable

    def connections_on(self, gid):
        return []

    def event_generators(self, gid):
        sched = arbor.explicit_schedule([1]) # one event at 1 ms
        return [arbor.event_generator('syn', 0, sched)]

    def probes(self, gid):
        return [arbor.cable_probe_point_state(target=0, mechanism='Ampa', state='g_active')]

    def global_properties(self, kind):
        return self.props

def run_sim(rec):
    context = arbor.context()
    decomp = arbor.partition_load_balance(rec, context)
    sim = arbor.simulation(rec, decomp, context)

    gid = 0
    sim.record(arbor.spike_recording.all)
    handles = sim.sample((gid, 0), arbor.regular_schedule(0.025))
    sim.run(5)
    samples, meta = sim.samples(handles)[0]

    return samples

def run_ground_truth():
    return run_sim(single_cell_recipe('gt'))

def run_approximation(a, b, c):
    return run_sim(single_cell_recipe('approx', a, b, c))

data_ap = run_approximation(0, 0.0, 0.00)
#data_gt = run_ground_truth()

import matplotlib.pyplot as plt
plt.plot(data_ap[:,0], data_ap[:, 1]) 
plt.xlabel('t(ms)')
plt.ylabel('Trelease(M)')
plt.title('approx plot')
plt.savefig('approx.png')

#plt.clf()
#plt.plot(data_gt[:,0], data_gt[:, 1]) 
#plt.xlabel('t(ms)')
#plt.ylabel('Trelease(M)')
#plt.title('ground truth plot')
#plt.savefig('gt.png')
