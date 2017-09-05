import numpy as np
from scipy.sparse import find as find_alive


class ConwayUniverseError(Exception):

    pass


class ConwayIndexError(ConwayUniverseError, IndexError):

    pass


class ConwayUniverse(object):

    ERR_WIDTH = "Given index ({}) exceeds the width of the Universe."
    ERR_HEIGHT = "Given index ({}) exceeds the height of the Universe."

    INFO_DEAD = ""

    def __init__(self, width, height=0, cycles=0, boundary=True,
                 quiet=False):

        self.__width = abs(int(width))

        if height == 0:
            self.__height = abs(int(width))
        else:
            self.__height = abs(int(height))

        self.__cycles = cycles

        self.__boundary = boundary

        self.__quiet = quiet

        # Let's create an empty universe.
        # There's no life...
        self.unimat = self.reshape()

    @property
    def width(self):
        """
        Defines the width of the Universe.
        """

        return self.__width

    @width.setter
    def width(self, width):
        self.__width = width

    @property
    def height(self):
        """
        Defines the height of the Universe.
        """

        return self.__height

    @height.setter
    def height(self, height):
        self.__height = height

    @property
    def cycles(self):
        """
        Number of cycles in the simulation.
        """

        return self.__cycles

    @cycles.setter
    def cycles(self, cycles):
        self.__cycles = cycles

    @property
    def boundary(self):
        """
        Determines if boundary conditions are used.
        Defaults to True.
        """

        return self.__boundary

    @boundary.setter
    def boundary(self, boundary):
        self.__boundary = boundary

    @property
    def quiet(self):
        """
        Determines if the class shows any notifications.
        """

        return self.__quiet

    @quiet.setter
    def quiet(self, quiet):
        self.__quiet = quiet

    def __getitem__(self, key):
        try:
            return self.unimat.__getitem__(key)
        except IndexError:

            if key[0] + 1 > self.width:
                message = self.ERR_WIDTH.format(key[0])
            elif key[1] + 1 > self.height:
                message = self.ERR_HEIGHT.format(key[1])

            raise IndexError(message)

    def __setitem__(self, key, alive):

        try:
            self.unimat.__setitem__(key, alive)

        except IndexError:

            if key[0] + 1 > self.width:
                message = self.ERR_WIDTH.format(key[0])
            elif key[1] + 1 > self.height:
                message = self.ERR_HEIGHT.format(key[1])

            raise ConwayIndexError(message)

    def __str__(self):
        return "Universe size is {} x {} cells.".format(
                self.width, self.height)

    def __translate_coordinate(self, coord, dimension):
        """
        """

        if not self.boundary:
            return coord

        if coord < 0:
            return dimension + coord
        elif coord > dimension - 1:
            return coord - dimension
        else:
            return coord

    def __middle(self, arg):
        """
        Returns index of the central element of given arg.
        """

        return int((len(arg) - 1) / 2)

    def reshape(self, width=0, height=0):

        """ Reshapes the Universe matrix."""

        if width == 0:
            w = self.width
        else:
            w = width

        if height == 0:
            h = self.height
        else:
            h = height

        return np.zeros(shape=(w, h), dtype=int)

    def nhood_moore(self, x, y, depth):

        """
        Returns Moore's neighbourhood of given cell.
        """

        actual_x = None
        actual_y = None

        moore_x = -1
        moore_y = -1

        dimension = 2 * depth + 1
        moore = self.reshape(dimension, dimension)

        for i in range(x - depth, x + depth + 1):

            moore_x += 1
            moore_y = -1

            actual_x = self.__translate_coordinate(i, self.width)

            for j in range(y - depth, y + depth + 1):

                moore_y += 1

                actual_y = self.__translate_coordinate(j, self.height)

                if self.boundary or (actual_x == i and actual_y == j):
                    moore[moore_x, moore_y] = self[actual_x, actual_y]
                else:
                        moore[moore_x, moore_y] = -1

        return moore

    def nhood_vN(self, x, y, depth):

        """
        Returns von Neuman's neighbourhood of given cell.
        """

        # Starting from Moore's:
        nh = self.nhood_moore(x, y, depth)

        # Find the cell in the middle.
        middle = self.__middle(nh)

        # TODO
        # Rewrite using slices???
        for i in range(0, len(nh)):
            for j in range(0, len(nh)):

                if i != middle and j != middle:
                    nh[i, j] = -1

    def nhood(self, x, y, depth=1, type="MOORE"):
        """
        Returns a matrix of cells neighbouring a cell at given x and y.
        """

        if type == "MOORE":
            return self.nhood_moore(x, y, depth)
        elif type == "vN":
            return self.nhood_vN(x, y, depth)

    def count_nhood(self, nhood=None, alive=True):

        """
        Counts either live or dead cells in a given neighbourhood.
        """

        # Count life in the whole Universe...
        if nhood is None:
            nhood = self[:]

        # or just in a given nhood.
        count = (nhood == int(alive)).sum()

        # Find out if the central element of the nhood is live or dead.
        # It's important because an alive central cell cannot be included
        # in its neighbourhood, thus count has to be decreased by 1.

        middle = self.__middle(nhood)

        if sum(len(nh) for nh in nhood) < 9 and nhood[middle, middle] == 1:
            count -= 1

        return count

    def seed(self):
        """
        Set up the initial state of the Universe.
        """

        pass

    def run(self):
        """
        Run the Universe.
        Uses find() method of scipy.sparse module to identify alive cells.

        Using nstate(x, y, alive = False), dead cells neighbouring every alive
        cell are going to be identified. These are going to be checked against
        their alive neighbourhood, so the last rule can be verified.
        """

        self.seed()

        for cycle in range(1, self.cycles + 1):

            if not self.quiet:
                print("Cycle #%i." % (cycle))

            if len(find_alive(self.unimat[:])[0]) == 0:
                if not self.quiet:
                    print("Universe is dead! Aborting simulation.")

                return

            ucopy = self.reshape()
            cells_alive = find_alive(self.unimat)

            for position in range(0, len(cells_alive[0])):

                x = cells_alive[0][position]
                y = cells_alive[1][position]

                nh = self.nhood(x, y)
                nalive = self.count_nhood(nh)
                middle = self.__middle(nh)

                if not self.quiet:
                    print("(%i, %i): %i " % (x, y, nalive))

                # Let's find out if the conditions for alive cells are met.
                actual_x = self.__translate_coordinate(x, self.width)
                actual_y = self.__translate_coordinate(y, self.height)
                ucopy[x, y] = int(nalive == 2 or nalive == 3)

                # Find all dead neighbours of current (alive) cell.
                for xd in range(0, len(nh)):
                    for yd in range(0, len(nh)):

                        if nh[xd, yd] == 0:
                            nh_dead = self.nhood(x, y)

                            if self.count_nhood(nh_dead) == 3:

                                actual_x = self.__translate_coordinate(
                                            x - middle + xd, self.width)

                                actual_y = self.__translate_coordinate(
                                            y - middle + yd, self.height)

                                ucopy[actual_x, actual_y] = True

                actual_x = self.__translate_coordinate(x, self.width)
                actual_y = self.__translate_coordinate(y, self.height)
                ucopy[actual_x, actual_y] = int(nalive == 3)

            self.unimat = ucopy
