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
        sched = arbor.explicit_schedule([1, 1.5, 2]) # one event at 1 ms
        return [arbor.event_generator('syn', 0, sched)]

    def probes(self, gid):
        if self.type == 'approx':
            #return [arbor.cable_probe_point_state(target=0, mechanism='Ampa', state='Trelease')]
            return [arbor.cable_probe_membrane_voltage('"synapse_site"')]
        elif self.type == 'gt':
            #return [arbor.cable_probe_point_state(target=0, mechanism='Ampa', state='Trelease')]
            return [arbor.cable_probe_membrane_voltage('"synapse_site"')]
        else: 
            return []

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

import matplotlib.pyplot as plt

#data_ap = run_approximation(3.01400226e-01, -2.06082296e-01, 1.80014261e-05)
#data_ap = run_approximation(3.55466060e-01, -2.16024506e-01,  1.19600159e-05)
#data_ap = run_approximation(3.91748572e-01, -2.19242639e-01,  9.90645981e-06)
#data_ap = run_approximation(4.24991430e-01, -2.09281013e-01,  8.30864218e-06)
#data_ap = run_approximation( 4.24842499e-01, -2.09230451e-01,  8.31083135e-06)
data_ap = run_approximation(4.81017695e-01, -3.56942029e-01, 7.81801390e-06)
data_gt = run_ground_truth()

plt.plot(data_ap[:,0], data_ap[:, 1], 'b') 
plt.plot(data_gt[:,0], data_gt[:, 1], 'g') 
plt.xlabel('t(ms)')
plt.ylabel('Trelease(M)')
plt.title('comparison plot')
plt.savefig('comp.png')
