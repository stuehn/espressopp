from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_List
class ListLocal(SetLocal, pairs_List):
    'Stores a list of pairs for a pairs::Computer to be applied.'
    def __init__(self, *args, **kywds):

        tuplen = len(args)
        dctlen = len(kywds)
        sumlen = tuplen + dctlen

        if sumlen == 2:
            if tuplen == 0:
                bc=kywds['bc']; storage=kywds['storage']
            elif tuplen == 1:
                bc=args[0]; storage=kywds['storage']
            else:
                bc=args[0]; storage=args[1]
            cxxinit(self, pairs_List, bc, storage)

        elif sumlen == 3:
            if tuplen == 0:
                bc=kywds['bc']; storage1=kywds['storage1']; storage2=kywds['storage2']
            elif tuplen == 1:
                bc=args[0]; storage1=kywds['storage1']; storage2=kywds['storage2']
            elif tuplen == 2:
                bc=args[0]; storage1=args[1]; storage2=kywds['storage2']
            else:
                bc=args[0]; storage1=args[1]; storage2=args[2]
            cxxinit(self, pairs_List, bc, storage1, storage2)
        else:
            raise ValueError('Number of arguments to constructor of pairs::List is invalid.')


if pmi.IS_CONTROLLER:
    class List(Set):
        'PMI class of a list of pairs'
        pmiproxydefs = \
            dict(cls = 'espresso.pairs.ListLocal', 
                 pmicall = [ 'addPair', 'deletePair', 'size', 'findPair' ])

        # TODO: addPair should check whether the ids are valid
