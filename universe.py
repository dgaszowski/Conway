"""
Implementation of Conway's Game Of Life.

Copyleft 2017
GNU GPL v. 3.0.

dgaszowski@gmail.com
"""

import numpy as np
from scipy.sparse import find as find_alive


class ConwayUniverseError(Exception):

    pass


class ConwayIndexError(ConwayUniverseError, IndexError):

    pass


class ConwayValueError(ConwayUniverseError, ValueError):

    pass


class ConwayUniverse(object):

    ERR_WIDTH = "Given index ({}) exceeds the width of the Universe."
    ERR_HEIGHT = "Given index ({}) exceeds the height of the Universe."
    ERR_DIMENSIONS = """Dimensions of the Universe have been already set.
    It is not allowed to change the dimensions during
    the simulation.
    """

    INFO_DEAD = ""

    CELL_ALIVE = 1
    CELL_DEAD = 0
    CELL_INACTIVE = -1

    def __init__(self, width, height=None, cycles=0, boundary=True,
                 quiet=False):

        self.__width = width

        if height is None:
            self.__height = width
        else:
            self.__height = height

        self.cycles = cycles

        self.boundary = boundary

        self.quiet = quiet

        # Let's create an empty universe an store it in a 2D matrix.
        self.__unimat = self.reshape()

    @property
    def width(self):
        """
        Returns the width of the Universe.
        """

        return self.__width

    @property
    def height(self):
        """
        Defines the height of the Universe.
        """

        return self.__height

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
        self.__boundary = bool(boundary)

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
            return self.__unimat.__getitem__(key)
        except IndexError:

            if not self.boundary:
                return

            if key[0] > self.width - 1:
                message = self.ERR_WIDTH.format(key[0])
            elif key[1] > self.height - 1:
                message = self.ERR_HEIGHT.format(key[1])

            raise ConwayIndexError(message)

    def __setitem__(self, key, alive):

        try:
            self.__unimat.__setitem__(key, alive)

        except IndexError:

            if key[0] > self.width - 1:
                message = self.ERR_WIDTH.format(key[0])
            elif key[1] > self.height - 1:
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

            if not self.boundary:
                if actual_x < 0 or actual_x > self.width - 1:
                    moore[moore_x, moore_y] = -1
                    continue

            for j in range(y - depth, y + depth + 1):

                moore_y += 1

                actual_y = self.__translate_coordinate(j, self.height)

                if not self.boundary:
                    if actual_y < 0 or actual_y > self.height - 1:
                        moore[moore_x, moore_y] = -1
                        continue

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

    def count_nhood(self, nhood=None):

        """
        Counts either live or dead cells in a given neighbourhood.
        """

        # Count life in the whole Universe...
        if nhood is None:
            nhood = self[:]

        # or just in a given nhood.

        middle_x = self.__middle(nhood[0])
        middle_y = self.__middle(nhood[1])
        count = (nhood == self.CELL_ALIVE).sum()

        if nhood[middle_x, middle_y] == 1:
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

            if len(find_alive(self.__unimat[:])[0]) == 0:
                if not self.quiet:
                    print("Universe is dead! Aborting simulation.")

                return

            ucopy = self.reshape()
            cells_alive = find_alive(self.__unimat)

            for position in range(len(cells_alive[0])):

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
                for xd in range(len(nh)):
                    for yd in range(len(nh)):

                        if nh[xd, yd] == 0:
                            nh_dead = self.nhood(x, y)

                            if self.count_nhood(nh_dead) == 3:

                                actual_x = self.__translate_coordinate(
                                            x - middle + xd, self.width)

                                actual_y = self.__translate_coordinate(
                                            y - middle + yd, self.height)

                                try:
                                    ucopy[actual_x, actual_y] = True
                                except IndexError:
                                    continue

                actual_x = self.__translate_coordinate(x, self.width)
                actual_y = self.__translate_coordinate(y, self.height)
                ucopy[actual_x, actual_y] = int(nalive == 3)

            self.__unimat = ucopy
