# when we're in pylab mode, the next two imports are not necessary
# we do it here for correctness sake, iow your code will also run without pylab mode

import numpy as np
import math
import matplotlib.pyplot as plt

import matplotlib.animation as animation

# coefficient of restitution (ratio of velocity after and before bounce)
# see http://en.wikipedia.org/wiki/Coefficient_of_restitution
#cor = 0.95
cor = 0.95

# bounds of the room
xlim = (0,1)
ylim = (0,1)

# coulomb factor
coulomb = 0.005

# electron and proton mass
m_e = 1
m_p = 1836.15267


delta_t = 0.002

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=xlim, ylim=ylim)
ax.grid()

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

class Ball():

    def __init__(self, m, q, xy, v0):
        """
        :param xy: Initial position.
        :param v: Initial velocity.
        """
        self.xy = np.array(xy)
        self.v = np.array(v0)
        self.m = m
        self.q = q
        
        print('markersize {}'.format( 5 * self.m**(1/3) ))
        color = 'red' if self.q < 0 else 'blue'

        self.scatter, = ax.plot([], [], color=color, marker='o', markersize = 15 * max(1, 0.2 * self.m**(1/3)))
        
    def get_force(self, balls):
        f = np.array((0.0, 0.0))
        for other in balls:
            if other is self:
                continue
            r = self.xy - other.xy
            dist2 = np.inner(r, r)
            f += self.q * other.q * coulomb * normalize(r) / dist2
        return f
        
    def first_verlet_step(self, balls):
        prev_xy = self.xy
        self.xy = prev_xy + delta_t * self.v + 0.5 * delta_t * delta_t * self.get_force( balls ) / self.m
        self.v = (self.xy - prev_xy) / delta_t
        #self.scatter.set_data(self.xy)
        
    def next_verlet_step(self, balls):
        prev_xy = self.xy
        self.xy = prev_xy + delta_t * self.v + delta_t * delta_t * self.get_force( balls ) / self.m
        self.v = (self.xy - prev_xy) / delta_t
        self.scatter.set_data(self.xy)

    def step(self, balls):
        if self.xy[0] <= xlim[0]:
            # hit the left wall, reflect x component
            self.v[0] = cor * np.abs(self.v[0])

        elif self.xy[0] >= xlim[1]:
            self.v[0] = - cor * np.abs(self.v[0])

        if self.xy[1] <= ylim[0]:
            # hit the bottom wall, reflect y component
            self.v[1] = cor * np.abs(self.v[1])

        elif self.xy[1] >= ylim[1]:
            self.v[1] = - cor * np.abs(self.v[1])
            
        self.next_verlet_step(balls)      

#    def update_position(self):
#        self.xy += delta_t * self.v
#
#        self.xy[0] = np.clip(self.xy[0], xlim[0], xlim[1])
#        self.xy[1] = np.clip(self.xy[1], ylim[0], ylim[1])
#        self.scatter.set_data(self.xy)
        
    def energy(self, balls):
        result = 0.5 * self.m * np.inner(self.v, self.v)
        for other in balls:
            if other is self:
                continue
            dist = np.linalg.norm(self.xy - other.xy)
            result += 0.5 * self.q * other.q * coulomb / dist
        return result


balls = [
    Ball(m_p,  1, (0.3,0.18), (0.023,0.012)),
    Ball(m_p,  1, (0.5,0.78), (-0.025,0.13)),
    Ball(m_e, -1, (0.04,0.67), (-0.045,0.043)),
    Ball(m_e, -1, (0.8,0.3), (-0.025,0.033))]
#total_energy = sum(ball.energy(balls) for ball in balls)

for ball in balls:
        ball.first_verlet_step(balls)

def init():
    return []

def animate(t):
    # t is time in seconds
    global xy, v
    
#    print('before 1 {}'.format(sum(ball.energy(balls) for ball in balls)))
    for ball in balls:
        ball.step(balls)
#    print('after 1   {}'.format(sum(ball.energy(balls) for ball in balls)))
    
    
    #print('before  {}'.format(sum(ball.energy(balls) for ball in balls)))
    #print('desired {}'.format(total_energy))
    #energies = [ball.energy(balls) for ball in balls]
    #fastest = max(balls, key=lambda b: np.inner(b.v, b.v))
    #prev_v = fastest.v
    #energy_change = total_energy - sum(energies)
    #factor = math.sqrt(max(0, 1 + energy_change / np.inner(prev_v, prev_v)))
#    print('energy fastest {}'.format(fastest.energy(balls)))
#    print('velocity fastest {}'.format(fastest.v))
#    print('energy_change {}'.format(energy_change))
#    print('factor {}'.format(factor))
    #fastest.v = factor * prev_v
#    print('after energy fastest {}'.format(fastest.energy(balls)))
#    print('after velocity fastest {}'.format(fastest.v))
    #print('after   {}'.format(sum(ball.energy(balls) for ball in balls)))
    
    # have to return an iterable
    return [ball.scatter for ball in balls]

# interval in milliseconds
# we're watching in slow motion (delta t is shorter than interval)
ani = animation.FuncAnimation(fig, animate, np.arange(0,100,delta_t), init_func=init, interval=10, blit=True)

plt.show()
