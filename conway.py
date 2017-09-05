# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from universe import ConwayUniverse as Universe

if __name__ == "__main__":

    u = Universe(16)
    print(u)

    u[13, 17] = 1
    u[12, 16] = 1
    u[13, 16] = 1
    u[0, 0] = 1

    u.cycles = 1000
    u.run()
    print("There are {} alive cells in the Universe.".format(u.count_nhood()))


