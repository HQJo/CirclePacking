from matplotlib.patches import Circle
import sys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(1)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    theta = np.linspace(0, 2*np.pi, 400)

    num = int(sys.argv[1])

    f = open('output.txt', 'r')
    data = f.read().split('\n')

    for d in data[0:num]:
        x0, y0, r0 = d.split(' ')
        x0, y0, r0 = float(x0), float(y0), float(r0)
        #ax.add_patch(Circle(xy=(x,y),radius=r, radius=0.1))
        x = x0 + r0 * np.cos(theta)
        y = y0 + r0 * np.sin(theta)
        ax.plot(x, y)
        #ax.plot(x0, y0)

    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)
    plt.show()
