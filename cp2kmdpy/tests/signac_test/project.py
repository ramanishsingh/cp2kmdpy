from flow import FlowProject, directives
import os

@FlowProject.label
def volume_computed(job):
    return job.isfile("volume.txt")


@FlowProject.operation
@FlowProject.post(volume_computed)
@directives(nranks=4)
def compute_volume(job):
    volume = job.sp.N * job.sp.kT / job.sp.p
    #for i in range(1000000000000):
    #    if i%100000000==0:
    #        print(i)
    with open(job.fn('volume.txt'), 'w') as file:
        file.write(str(volume) + '\n')

@FlowProject.operation
def printing(job):
    print(volume_computed)
if __name__ == '__main__':
    FlowProject().main()
