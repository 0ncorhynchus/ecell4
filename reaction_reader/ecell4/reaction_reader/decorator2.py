import types
import numbers
import copy
import functools

import species
import decorator
import parseobj


def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        sp = species.Species()
        for elem in obj._elements():
            su = species.Subunit(elem.name)
            if elem.args is not None:
                for mod in elem.args:
                    if (not (isinstance(mod, parseobj.ParseObj)
                            or isinstance(mod, parseobj.AnyCallable))
                        or mod._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (mod))
                    arg = mod._elements()[0]
                    name, binding = arg.name, arg.modification
                    if binding is None:
                        su.add_modification(name, "", "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, "", binding)
            if elem.kwargs is not None:
                for name, value in elem.kwargs.items():
                    if (not (isinstance(value, parseobj.ParseObj)
                            or isinstance(value, parseobj.AnyCallable))
                        or value._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (value))
                    arg = value._elements()[0]
                    state, binding = str(arg.name), arg.modification
                    if binding is None:
                        su.add_modification(name, state, "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, state, binding)
            sp.add_subunit(su)
        return (sp, )
    elif isinstance(obj, parseobj.InvExp):
        return (None, )
    elif isinstance(obj, parseobj.AddExp):
        subobjs = obj._elements()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)

    raise RuntimeError, 'invalid expression; "%s" given' % str(obj)

def generate_ReactionRule(lhs, rhs, opts):
    # if len(lhs) == 0:
    #     if len(rhs) != 1:
    #         raise RuntimeError, (
    #             "the number of products must be 1; %d given" % len(rhs))
    #     return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    if len(lhs) == 0 or len(lhs) == 1 or len(lhs) == 2:
        return species.ReactionRule(lhs, rhs, opts)
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))

def generate_Option(opt):
    # if not (isinstance(opt, parseobj.AnyCallable)
    #     or isinstance(opt, parseobj.ParseObj)):
    #     raise RuntimeError

    if opt._size() != 1:
        raise RuntimeError

    elem = opt._elements()[0]
    if elem.name == "IncludeReactants" or elem.name == "ExcludeReactants":
        if not (len(elem.args) == 2
            and type(elem.args[0]) == int
            and (isinstance(elem.args[1], parseobj.AnyCallable)
                or isinstance(elem.args[1], parseobj.ParseObj))):
            raise RuntimeError

        if isinstance(elem.args[1], parseobj.ParseObj):
            raise RuntimeError, "only a subunit name is allowed here."

        pttrn = elem.args[1]._elements()[0].name
        if elem.name == "ExcludeReactants":
            return (species.ExcludeReactants(elem.args[0], pttrn),
                species.ExcludeProducts(elem.args[0], pttrn))
        elif elem.name == "IncludeReactants":
            return (species.IncludeReactants(elem.args[0], pttrn),
                species.IncludeProducts(elem.args[0], pttrn))
    elif elem.name == "IncludeProducts" or elem.name == "ExcludeProducts":
        if not (len(elem.args) == 2
            and type(elem.args[0]) == int
            and (isinstance(elem.args[1], parseobj.AnyCallable)
                or isinstance(elem.args[1], parseobj.ParseObj))):
            raise RuntimeError

        if isinstance(elem.args[1], parseobj.ParseObj):
            raise RuntimeError, "only a subunit name is allowed here."

        pttrn = elem.args[1]._elements()[0].name
        if elem.name == "ExcludeProducts":
            return (species.ExcludeProducts(elem.args[0], pttrn),
                species.ExcludeReactants(elem.args[0], pttrn))
        elif elem.name == "IncludeProducts":
            return (species.IncludeProducts(elem.args[0], pttrn),
                species.IncludeReactants(elem.args[0], pttrn))
    else:
        # raise RuntimeError
        return (opt, None)
    return (opt, opt)

def generate_Options1(opts):
    retval = []
    for opt in opts:
        if (isinstance(opt, parseobj.AnyCallable)
            or isinstance(opt, parseobj.ParseObj)):
            lhs, rhs = generate_Option(opt)
            if lhs is not None:
                retval.append(lhs)
        # elif isinstance(opt, numbers.Number):
        #     retval.append(opt)
        # else:
        #     raise RuntimeError, "an invalid option [%s] given." % (opt)
        retval.append(opt)
    return retval

def generate_Options2(opts):
    retval1, retval2 = [], []
    for opt in opts:
        if (isinstance(opt, parseobj.AnyCallable)
            or isinstance(opt, parseobj.ParseObj)):
            lhs, rhs = generate_Option(opt)
            if lhs is not None:
                retval1.append(lhs)
            if rhs is not None:
                retval2.append(rhs)
        elif ((isinstance(opt, types.ListType)
            or isinstance(opt, types.TupleType))
            and len(opt) == 2):
            # if (isinstance(opt[0], numbers.Number)
            #     and isinstance(opt[1], numbers.Number)):
            #     raise RuntimeError
            retval1.append(opt[0])
            retval2.append(opt[1])
        else:
            raise RuntimeError, "an invalid option [%s] given." % (opt)
    return retval1, retval2

class SpeciesAttributesCallback(decorator.Callback):

    def __init__(self, *args):
        decorator.Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        if not isinstance(obj, parseobj.OrExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))
        elif len(obj._elements()) != 2:
            raise RuntimeError, 'only one attribute is allowed. [%d] given' % (
                len(obj._elements()))

        lhs, rhs = obj._elements()

        species_list = generate_Species(lhs)
        if len(species_list) != 1:
            raise RuntimeError, (
                'only a single species must be given; %d given'
                % len(species_list))

        sp = species_list[0]
        if sp is None:
            raise RuntimeError, "no species given [%s]" % (repr(obj))

        self.bitwise_operations.append((sp, rhs))

    def notify_comparisons(self, obj):
        raise RuntimeError, (
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(decorator.Callback):

    def __init__(self):
        decorator.Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_comparisons(self, obj):
        if not isinstance(obj, parseobj.CmpExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))
        elif isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(lhs, parseobj.OrExp):
            lhs = lhs._elements()[0]

        if isinstance(rhs, parseobj.OrExp):
            opts = rhs._elements()[1: ]
            rhs = rhs._elements()[0]
        else:
            opts = []

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp in lhs if sp is not None)
        rhs = tuple(sp for sp in rhs if sp is not None)

        if isinstance(obj, parseobj.EqExp) or isinstance(obj, parseobj.NeExp):
            opts1, opts2 = generate_Options2(opts)
            self.comparisons.append(generate_ReactionRule(lhs, rhs, opts1))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, opts2))
        elif isinstance(obj, parseobj.GtExp):
            opts = generate_Options1(opts)
            self.comparisons.append(generate_ReactionRule(lhs, rhs, opts))
        else:
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))

species_attributes = functools.partial(decorator.parse_decorator, SpeciesAttributesCallback)
reaction_rules = functools.partial(decorator.parse_decorator, ReactionRulesCallback)

class AnyCallableGenerator(dict):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

        self.__cache = decorator.Callback()

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        retval = self.get(key)
        if retval is None:
            return parseobj.AnyCallable(self.__cache, key)
        return retval

def create_species(expr):
    vardict = AnyCallableGenerator()
    obj = eval(expr, globals(), vardict)
    retval = generate_Species(obj)
    if len(retval) != 1:
        raise RuntimeError, "multiple species were defined in the expression."
    return retval[0]


if __name__ == "__main__":
    def cmp_subunit(sp, idx1, idx2):
        retval = species.CmpSubunit(sp)(sp.subunits[idx1], sp.subunits[idx2])
        print "%s %s %s in %s" % (
            str(sp.subunits[idx1]),
            "==" if retval == 0 else (">" if retval == 1 else "<"),
            str(sp.subunits[idx2]),
            str(sp))
        return retval

    sp1 = create_species("Grb2(SH2^1,SH3^2).Grb2(SH2^3,SH3^4).Grb2(SH2^5,SH3^6).Grb2(SH2^7,SH3^8).Shc(PTB^9,Y317=pY^3).Shc(PTB^10,Y317=pY^7).Sos(dom^2).Sos(dom^4).Sos(dom^6).Sos(dom^8).egf(r^11).egf(r^12).egfr(l^11,r^13,Y1068=pY^1,Y1148=pY^9).egfr(l^12,r^13,Y1068=pY^5,Y1148=pY^10)")
    sp2 = create_species("egf(r^1).egfr(l^1,r^4,Y1068=pY^2,Y1148=pY^6).Grb2(SH2^2,SH3).egf(r^3).egfr(l^3,r^4,Y1068=pY,Y1148=pY^9).Shc(PTB^6,Y317=pY^7).Grb2(SH2^7,SH3).Shc(PTB^9,Y317=pY^10).Grb2(SH2^10,SH3)")
    sp3 = create_species("A(bs^1).B(l^1,r^2).A(bs^2)")
    sp4 = create_species("A(bs^1).B(l^1,r^3).B(r^3,l^2).A(bs^2)")
    sp5 = create_species("A(bs^1).B(l^2,r^3).B(r^3,l^1).A(bs^2)")

    sp6 = create_species("L(l1^1,l2^2).L(l1^3,l2^4).L(l1^5,l2^6).R(r1^3,r2^2).R(r1^5,r2^4).R(r1^1,r2^6)")

    cmp_subunit(sp1, 0, 1)
    cmp_subunit(sp1, 0, 2)
    cmp_subunit(sp1, 0, 3)
    cmp_subunit(sp1, 10, 11)
    cmp_subunit(sp2, 0, 3)
    cmp_subunit(sp2, 2, 6)
    cmp_subunit(sp2, 2, 8)
    cmp_subunit(sp2, 6, 8)
    cmp_subunit(sp3, 0, 2)
    cmp_subunit(sp4, 0, 3)
    cmp_subunit(sp4, 1, 2)

    cmp_subunit(sp6, 0, 1)
    cmp_subunit(sp6, 0, 2)
    cmp_subunit(sp6, 1, 2)
    cmp_subunit(sp6, 3, 4)
    cmp_subunit(sp6, 3, 5)
    cmp_subunit(sp6, 4, 5)

    print ""
    sp1.sort()
    print sp1
    sp2.sort()
    print sp2
    sp3.sort()
    print sp3
    sp4.sort()
    print sp4
    sp5.sort()
    print sp5
    sp6.sort()
    print sp6
    print ""

    import random
    sp = create_species("L(l1^1,l2^2).L(l1^3,l2^4).L(l1^5,l2^6).R(r1^3,r2^2).R(r1^5,r2^4).R(r1^1,r2^6)")
    # sp = create_species("L(l1^1,l2^2).L(l1^3,l2^4).L(l1^5,l2^6).R(r1^1,r2^6).R(r1^3,r2^2).R(r1^5,r2^4)")
    # sp = create_species("L(l^1,r^2).L(l^2,r^3).L(l^3,r^1)")
    # sp = create_species("L(l^1,r^2).L(l^2,r^3).L(l^3,r^4).L(l^4,r^1)")
    newbs = range(1, 7)
    # print 'ORIGINAL    :', sp
    for _ in range(10):
        random.shuffle(sp.subunits)
        random.shuffle(newbs)
        sp.update_indices()
        for su in sp.subunits:
            for mod in su.modifications.keys():
                state, bs = su.modifications[mod]
                if bs.isdigit():
                    su.modifications[mod] = (state, str(newbs[int(bs) - 1]))
        # print '[%d] SHUFFLED:' % _, sp
        sp.sort()
        print '[%d] SORTED  :' % _, sp

    print ""
    sp1 = create_species("L(l1^3,l2^5).R(r1^6,r2^5).L(l1^6,l2^2).L(l1^4,l2^1).R(r1^3,r2^1).R(r1^4,r2^2)")
    sp1.sort()
    print sp1
    sp2 = create_species("L(l1^3,l2^5).R(r1^1,r2^6).R(r1^4,r2^5).R(r1^3,r2^2).L(l1^1,l2^2).L(l1^4,l2^6)")
    sp2.sort()
    print sp2

    print ""
    sp1 = create_species("L(l^2,r^3).L(l^3,r^1).L(l^1,r^2)")
    print sp1
    sp1.sort()
    print sp1
    sp2 = create_species("L(l^1,r^3).L(l^2,r^1).L(l^3,r^2)")
    print sp2
    sp2.sort()
    print sp2
