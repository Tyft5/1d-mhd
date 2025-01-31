import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
from sys import argv
import re

fnames = argv[1::5]

tl_arr = [int(re.search('\d+', nm)[0]) for nm in fnames]

fnames = [f for _, f in sorted(zip(tl_arr, fnames))]

tl_arr = sorted(tl_arr)

tl_meshes = []
for fn in fnames:
    mesh = []
    for line in open(fn):
        mesh.append([float(n) for n in line.split(',')])
    tl_meshes.append(mesh)

tl_meshes = np.array(tl_meshes)
x_mesh = np.linspace(0, 1, 800)

animate = len(tl_meshes) > 1

outfn = fn.strip('.out') + '.png'
mesh = tl_meshes[0]

fig, axes = plt.subplots(2, 3, sharex='all')
#fig.set_size_inches(10, 8)
fig.set_dpi(150)

toplot = mesh[:, [0, 1, 2, 4, 6, 7]]
line_arr = list(range(6))

ylims = [
    (0.1, 1.1),
    (-0.3, 0.7),
    (0, 1.6),
    (-1, 1),
    (0.1, 1),
    (-1, 0.5)
]
titles =[ r"$\rho$", r"$u_x$", r"$u_y$", r"$B_y$", r"$e$", r"$p$" ]

for i, ax in enumerate(axes.reshape(-1)):
    line_arr[i] = ax.scatter(x_mesh, toplot[:, i], s=0.5)
    ax.set_ylim(ylims[i])
    ax.set_ylabel(titles[i])

plt.tight_layout()
# if not animate:
    # plt.suptitle(f"TL {tl_arr[0]}")
    # plt.subplots_adjust(top=0.9)

def update(frame):
    mesh = tl_meshes[frame]
    toplot = mesh[:, [0, 1, 2, 4, 6, 7]]
    for i in range(6):
        line_arr[i].set_offsets(np.array([x_mesh, toplot[:, i]]).T)
    
    plt.suptitle(f"TL {tl_arr[frame]}")
    
    return (line_arr,)

if animate:
    ani = anim.FuncAnimation(fig=fig, func=update, frames=len(tl_meshes), interval=200)

#    ani.save('brio_wu_anim.mp4', fps=5, extra_args=['-vcodec', 'libx264'])
else:
    plt.savefig(outfn, bbox_inches='tight', dpi=200)
plt.show()
