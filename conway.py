# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from universe import ConwayUniverse as Universe

if __name__ == "__main__":

    u = Universe(1000)
    u.boundary = True
    print(u)

    u[13, 17] = 1
    u[12, 16] = 1
    u[1, 2] = 1
    u[1, 1] = 1
    u[1, 0] = 1
    u[0, 0] = 1

    u.cycles = 10000
    u.run()
    print("There are {} alive cells in the Universe.".format(u.count_nhood()))


