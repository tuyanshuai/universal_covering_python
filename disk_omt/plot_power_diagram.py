"""
function plot_power_diagram(pd)
% plot cells of power diagram, can be slow
gcf;
hold on
bd = compute_bd(pd.face);
disk = pd.uv(bd,:);
% plot(disk(:,1),disk(:,2),'r-')
for i = 1:length(pd.cell)
    pi = pd.dpe(pd.cell{i},:);
    plot(pi(:,1),pi(:,2),'b-','Linewidth',1.5);
end
box = [min(pd.uv(:,1)),max(pd.uv(:,1)),min(pd.uv(:,2)),max(pd.uv(:,2))];
axis equal
axis(box)
"""


"""
Plot_mesh
"""
import numpy as np
import matplotlib.pyplot as plt
from  algebra import *

# from algebra import *
def plot_power_diagram(pd):
    bd = compute_bd(pd["face"])
    bd = np.append(bd, bd[0])
    disk = pd["uv"][bd,:]

    fig = plt.figure()
    plt.plot(disk[:, 0], disk[:, 1], 'r-')
    plt.plot(pd["uv"][:,0], pd["uv"][:,1],'g.')




    for i in range(len(pd["cell"])):
        pi = pd["dpe"][pd["cell"][i].flatten(),:]
        plt.plot(pi[:, 0], pi[:, 1], 'b-', linewidth = 1.5)

    plt.axis('equal')
    plt.show()
