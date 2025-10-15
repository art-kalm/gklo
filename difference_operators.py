r"""
Difference Operators

AUTHORS:

- Artem Kalmykov (2025-10-14): Initial version
"""

# ****************************************************************************
#  Copyright (C) 2025 Artem Kalmykov <artem.o.kalmykov at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.sets.family import Family


class DifferenceOperators(CombinatorialFreeModule):
    r"""
    Algebra of localized difference operators `A`.

    Let `R` be a ring and `hbar` be its element. Let 
    
    .. MATH::
        F = R(x_1, ..., x_l, p_{i_1}, ..., p_{i_k}) 
    
    be the ring of rational functions in *functional variables* `p_{i_1}`, ..., `p_{i_k}` and *additional variables* `x_1`, ..., `x_l`; 
    here, `I = {i_1, ..., i_k}` is some index set and `l` is some number. Then, as a left F-module, 
    the algebra of (localized) difference operators is 
    
    .. MATH:: 
        A = F[u_{i_1}^{\pm 1}, ..., u_{i_k}^{\pm 1}], 
        
    with relations

    .. MATH::
        u_{i_a}^{\pm} f(x_1, ..., x_k, p_{i_1}, ... p_{i_a}, ..., p_{i_k}) = f(x_1, ..., x_k, p_{i_1}, ... p_{i_a} \pm hbar, ..., p_{i_k}) u_{i_a}^{\pm}, 

    and `u_{i_a}^{\pm}` commute among themselves; here, `f(x_1, ..., x_l, p_{i_1}, ..., p_{i_k})` is any function in `F`.

    INPUT::

    - ``base_ring`` -- the base ring `R`
    - ``I`` -- the index set
    - ``shift`` -- an element of ``base_ring``, corresponds to `hbar`
    - ``variable_prefix`` -- (default: ``p``) the prefix of the functional variables. The latter are represented as ``variable_prefix` + `_i`` where `i=str(i)` for `i` in `I` 
        (see commentary below)
    - ``differece_operator_prefix`` -- (default: ``u``) the prefix of difference operators. The latter are represented as ``difference_operator_prefix` + `_i`` where `i=str(i)` for `i` in `I` 
        (see commentary below)
    - ``additional_variables`` -- (default: ``()``) additional variables to which we do not assign any difference operators.

    .. NOTE:: if an element `i` of the index set ``I`` is iterable (for instance, of type `list`, `tuple`, `set`, etc), then the text representation
        of the corresponding functional variable is ```variable_prefix` + `_a[0]` + `_a[1]` + ... + `_a[m]`, where `a[i]` are values of iterator of `i`. 
        See examples below.

    .. NOTE:: to simplify the code (in fact, significantly: otherwise, event simplest functions like ``_add_()`` or an obvious coercion
        from `F` to `A` need to be redefined), the function ``self.base_ring()`` returns ``F`` instead of ``base``,
        although `A` is *not* an algebra over `F`, just a left module. 
        Also, for coercion purposes, it is convenient to set the category of ``self`` as ``AlgebrasWithBasis(F)`` which is, strictly speaking, not entirely correct.
        In my tests, it did not lead to any bugs, however, please
        let me know if it needs to be changed.

    EXAMPLES::

        sage: I = range(4); D = DifferenceOperators(QQ, I); D
        Localized algebra of difference operators over Rational Field with shift parameter 1 in functional variables (p_0, p_1, p_2, p_3), difference operators (u_0, u_1, u_2, u_3), and spectral variables ()
        sage: p = D.functional_variables(); p
        Finite family {0: p_0, 1: p_1, 2: p_2, 3: p_3}
        sage: u = D.difference_operators(); u
        Finite family {0: u_0, 1: u_1, 2: u_2, 3: u_3}
        sage: p[0]*u[0]
        p_0*u_0
        sage: u[0]*p[0]
        (p_0+1)*u_0
        sage: u[0]**(-1)*p[0]
        (p_0-1)*u_0^-1
        sage: u[0]*p[2]
        p_2*u_0
        sage: (p[0]*u[0])**(-1)
        (1/(p_0-1))*u_0^-1
        sage: a = p[0]*u[0]; a**(-1)
        (1/(p_0-1))*u_0^-1
        sage: a*a**(-1)
        1

        sage: P = PolynomialRing(QQ,'h'); h = P.gen()
        sage: I = [(i,j) for i in range(3) for j in range(2)]; DiffOps = DifferenceOperators(base_ring=P, I = I, shift = 2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s')); DiffOps
        Localized algebra of difference operators over Univariate Polynomial Ring in h over Rational Field with shift parameter 2*h in functional variables (w_0_0, w_0_1, w_1_0, w_1_1, w_2_0, w_2_1), difference operators (D_0_0, D_0_1, D_1_0, D_1_1, D_2_0, D_2_1), and spectral variables ('t', 's')
        sage: w = DiffOps.functional_variables(); w
        Finite family {(0, 0): w_0_0,  (0, 1): w_0_1,  (1, 0): w_1_0,  (1, 1): w_1_1,  (2, 0): w_2_0,  (2, 1): w_2_1}
        sage: d = DiffOps.difference_operators(); d
        Finite family {(0, 0): D_0_0,  (0, 1): D_0_1,  (1, 0): D_1_0,  (1, 1): D_1_1,  (2, 0): D_2_0,  (2, 1): D_2_1}
        sage: w[0,1]*d[0,1]
        w_0_1*D_0_1
        sage: d[0,1]*w[0,1]
        (w_0_1+2*h)*D_0_1
        sage: d[0,1]^(-1)*w[0,1]
        (w_0_1+(-2*h))*D_0_1^-1
        sage: DiffOps.additional_variables()
        (t, s)
        sage: s = DiffOps.additional_variables()[1]
        sage: d[2,0]*s
        s*D_2_0
        sage: s*d[2,0] == d[2,0]*s
        True
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, I, shift=1, variable_prefix='p', difference_operator_prefix='u', additional_variables=()):
        r"""
        Return the correct parent based upon input.
        Preprocesses the types of ``I`` and ``additional_variables`` to `tuple` to ensure a unique and hashable representation.

        EXAMPLES::

        sage: I = range(4); D = DifferenceOperators(QQ, I); D
        Localized algebra of difference operators over Rational Field with shift parameter 1 in functional variables (p_0, p_1, p_2, p_3), difference operators (u_0, u_1, u_2, u_3), and spectral variables ()
        sage: D1 = DifferenceOperators(QQ, [0,1,2,3])
        sage: D1 is D
        True
        sage: D2 = DifferenceOperators(QQ, [0,1,2,3], additional_variables='t')
        sage: D2 is D
        False
        """
        I = MonadTuple(I)
        additional_variables = MonadTuple(additional_variables)
        return super().__classcall__(cls, base_ring, I, shift, variable_prefix, difference_operator_prefix, additional_variables)

    def __init__(self, base_ring, I, shift, variable_prefix, difference_operator_prefix, additional_variables):
        self._shift = shift
        self._I = I
        self.difference_operator_prefix = difference_operator_prefix
        self._spectral_variables_str = additional_variables
        self.variable_prefix = variable_prefix
        diff_op_variables = [variable_prefix + '_' + _to_str(i) for i in I]
        all_variables = list(additional_variables) + diff_op_variables
        functional_base = PolynomialRing(base_ring, all_variables).fraction_field()

        basis_keys = IndexedFreeAbelianGroup(indices=I,bracket=False)
        category = AlgebrasWithBasis(functional_base)
        CombinatorialFreeModule.__init__(self, functional_base, basis_keys, prefix=difference_operator_prefix, category=category)

    def _left_module_base_ring(self):
        r"""
        Return `F` from the description, that is, the base of difference operators as a left module.

        .. NOTE:: so far, equivalent to ``self.base_ring()``; however, since the latter is technically not the base ring,
            I prefer to keep for future possible changes.

        EXAMPLES::

        sage: I = range(4); D = DifferenceOperators(QQ, I)
        sage: D._left_module_base_ring()
        Fraction Field of Multivariate Polynomial Ring in p_0, p_1, p_2, p_3 over Rational Field
        sage: D1 = DifferenceOperators(QQ, I, additional_variables='t')
        sage: D1._left_module_base_ring()
        Fraction Field of Multivariate Polynomial Ring in t, p_0, p_1, p_2, p_3 over Rational Field
        """
        return super().base_ring()

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

        sage: I = range(4); D = DifferenceOperators(QQ, I); D
        Localized algebra of difference operators over Rational Field with shift parameter 1 in functional variables (p_0, p_1, p_2, p_3), difference operators (u_0, u_1, u_2, u_3), and spectral variables ()
        """
        return "Localized algebra of difference operators over {} with shift parameter {} in functional variables {}, difference operators {}, and spectral variables {}".format(self.base_ring(), self._shift, tuple(self.functional_variables()), tuple(self.difference_operators()), self._spectral_variables_str)

    def shift(self):
        r"""
        Return ``shift`` parameter of ``self``.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I); D
            Localized algebra of difference operators over Rational Field with shift parameter 1 in functional variables (p_0, p_1, p_2, p_3), difference operators (u_0, u_1, u_2, u_3), and spectral variables ()
            sage: D.shift()
            1

            sage: I = [(i,j) for i in range(3) for j in range(2)]
            sage: DiffOps = DifferenceOperators(base_ring=P, I = I, shift = 2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s')); DiffOps
            Localized algebra of difference operators over Univariate Polynomial Ring in h over Rational Field with shift parameter 2*h in functional variables (w_0_0, w_0_1, w_1_0, w_1_1, w_2_0, w_2_1), difference operators (D_0_0, D_0_1, D_1_0, D_1_1, D_2_0, D_2_1), and spectral variables ('t', 's')
            sage: DiffOps.shift()
            2*h
        """
        return self._shift
        
    @cached_method
    def one(self):
        r"""
        Return 1 of ``self``.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ,I)
            sage: D.one()
            1
        """
        return self._from_dict({self._indices.one(): self._left_module_base_ring().one()})

    @cached_method
    def additional_variables(self):
        r"""
        Return additional variables of ``self`` to which we do *not* assign difference operators.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: D.additional_variables()
            (t,)
            sage: t = D.additional_variables()[0]; u = D.difference_operators()
            sage: t*u[0]
            t*u_0
            sage: u[0]*t
            t*u_0

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables=('t','s'))
            sage: D.additional_variables()
            (t, s)
        """
        d = self._left_module_base_ring().gens_dict()
        return tuple(d[i] for i in self._spectral_variables_str)

    @cached_method
    def functional_variables(self):
        r"""
        Return functional variables of ``self`` parameterized by ``I`` to which we assign difference operators.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: p = D.functional_variables(); p
            Finite family {0: p_0, 1: p_1, 2: p_2, 3: p_3}
            sage: p[0]**(-3)*p[1]^2*p[2]
            p_1^2*p_2/p_0^3
        """
        d = self._left_module_base_ring().gens_dict()
        return Family({i: d[self.variable_prefix + '_' + _to_str(i)] for i in self._I},
                      name="generator")

    @cached_method
    def difference_operators(self):
        r"""
        Return difference operators of ``self`` parameterized by ``I``.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: u = D.difference_operators(); u
            Finite family {0: u_0, 1: u_1, 2: u_2, 3: u_3}
            sage: p = D.functional_variables(); p
            Finite family {0: p_0, 1: p_1, 2: p_2, 3: p_3}
            sage: u[0]^(-1)*p[0]
            (p_0-1)*u_0^-1
            sage: u[1]^2*p[1]
            (p_1+2)*u_1^2
            sage: u[2]*p[3] == p[3]*u[2]
            True

            sage: D1 = DifferenceOperators(P, I, shift=-2*h, additional_variables='t')
            sage: u = D1.difference_operators(); u
            Finite family {0: u_0, 1: u_1, 2: u_2, 3: u_3}
            sage: p = D1.functional_variables()
            sage: u[0]^(-1)*p[0]
            (p_0+(2*h))*u_0^-1

            sage: I2 = [(i,j) for i in range(3) for j in range(2)]; I2
            [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
            sage: D2 = DifferenceOperators(base_ring=P, I=I2, shift=2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s'))
            sage: D2.difference_operators()
            Finite family {(0, 0): D_0_0,  (0, 1): D_0_1,  (1, 0): D_1_0,  (1, 1): D_1_1,  (2, 0): D_2_0,  (2, 1): D_2_1}
        """
        return Family({i : self.gen(i) for i in self._I})
    
    def gen(self, x):
        r"""
        Return the difference operator corresponding to ``x``.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: D.gen(0)
            u_0

            sage: I2 = [(i,j) for i in range(3) for j in range(2)]; I2
            [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
            sage: D2 = DifferenceOperators(base_ring=P, I=I2, shift=2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s'))
            sage: u = D2.gen((2,0)); u
            D_2_0
            sage: p = D2.functional_variables()[2,0]
            sage: u*p
            (w_2_0+2*h)*D_2_0
        """
        m = self._indices.gen(x)
        return self.element_class(self, {m: self.base_ring().one()})
    
    def _repr_term(self, m):
        r"""
        Return a string representation of the basis element (of `A` considered as a left module over `F`) indexed by ``m``.

        EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: D._repr_term(D.indices().gen(0)**2*D.indices().gen(1)**(-1))
            'u_0^2*u_1^-1'
        """
        if m == self.indices().one():
            return '1'
        return '*'.join(self.difference_operator_prefix + '_' + _to_str(i)
                        + ('^{}'.format(exp) if exp != 1 else '')
                        for i,exp in m._sorted_items())
    

    class Element(CombinatorialFreeModule.Element):
        def __call__(self, **kwargs):
            r"""
            Return ``self`` with substituted *additional* variables according to ``kwargs``.

            EXAMPLES::

                sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables=['t','s'])
                sage: t, s = D.additional_variables()
                sage: u = D.difference_operators(); p = D.functional_variables()
                sage: element = u[0]*(t-p[0]); element
                (t-p_0-1)*u_0
                sage: element(t=2*s)
                (2*s-p_0-1)*u_0
                sage: element(t=p[0])
                -u_0
            """
            if not kwargs.keys() <= set(self.parent()._spectral_variables_str):
                raise TypeError("Evaluation is defined only for additional variables")
            new_element_dict = {}
            for monom,coeff in self.monomial_coefficients().items():
                new_element_dict[monom] = coeff(**kwargs)
            return self.parent()._from_dict(new_element_dict)
        
        def _pow_int(self, deg):
            r"""
            Return ``deg`` power of ``self``. Defined *only* for elements having single monomial.

            EXAMPLES::

            sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
            sage: element = p[0]*u[0]; element
            p_0*u_0
            sage: element**2
            (p_0^2+p_0)*u_0^2
            sage: element^(-2)
            (1/(p_0^2-3*p_0+2))*u_0^-2
            sage: element**2*element^(-2)
            1
            """
            parent = self.parent()
            if deg == 0:
                return parent.one()
            if deg > 1:
                new_element = parent.one()
                for _ in range(deg):
                    new_element = new_element*self
                return new_element
            if deg < 0:
                if len(self.monomial_coefficients().items()) != 1:
                    raise TypeError("Element should be monomial")
                monom = next(iter(self.monomial_coefficients().keys()))
                aux_element = parent.monomial(monom**(-1))
                inverse = aux_element*self.monomial_coefficients()[monom]**(-1)
                new_element = parent.one()
                for _ in range(-deg):
                    new_element = new_element*inverse
                return new_element

        def _rmul_(self, other):
            r"""
            Return ``other``*``self`` where ``other`` is an element of `F`.

            EXAMPLES::

                sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
                sage: u = D.difference_operators(); p = D.functional_variables()
                sage: p[0]*u[0]
                p_0*u_0
                sage: p[0]*u[2]
                p_0*u_2

                sage: I2 = [(i,j) for i in range(3) for j in range(2)]; I2
                [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
                sage: D2 = DifferenceOperators(base_ring=P, I=I2, shift=2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s'))
                sage: t, s = D2.additional_variables()
                sage: w = D2.functional_variables(); d = D2.difference_operators()
                sage: t*w[1,0]*d[0,1]
                t*w_1_0*D_0_1
            """
            return super()._rmul_(other)
        
        def _lmul_(self,other):
            r"""
            Return ``self``*``other`` where ``other`` is an element of `F`.

            EXAMPLES::

                sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
                sage: u = D.difference_operators(); p = D.functional_variables()
                sage: u[0]*p[0]
                (p_0+1)*u_0
                sage: u[0]*p[2]
                p_2*u_0

                sage: I2 = [(i,j) for i in range(3) for j in range(2)]; I2
                [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
                sage: D2 = DifferenceOperators(base_ring=P, I=I2, shift=2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s'))
                sage: t, s = D2.additional_variables()
                sage: w = D2.functional_variables(); d = D2.difference_operators()
                sage: d[0,1]*(t*w[1,0])
                t*w_1_0*D_0_1
            """
            parent = self.parent()
            new_element_dict = {}
            for monom,coeff in self.monomial_coefficients().items():
                shifts = {parent.variable_prefix + '_' + _to_str(index): parent.functional_variables()[index] + deg*parent._shift for index,deg in monom.dict().items()}
                new_element_dict[monom] = coeff*other.subs(**shifts)
            return parent._from_dict(new_element_dict)
        
        def _mul_(self,other):
            r"""
            Return ``self``*``other`` where ``other`` is an element of `A`.

            EXAMPLES::

                sage: I = range(4); D = DifferenceOperators(QQ, I, additional_variables='t')
                sage: u = D.difference_operators(); p = D.functional_variables()
                sage: element1 = p[0]*u[0]; element1
                p_0*u_0
                sage: element1 = p[0]*u[1]; element1
                p_0*u_1
                sage: element2 = p[1]*u[2]; element2
                p_1*u_2
                sage: element1*element2
                (p_0*p_1+p_0)*u_1*u_2
                sage: element2*element1
                p_0*p_1*u_1*u_2

                sage: I2 = [(i,j) for i in range(3) for j in range(2)]; I2
                [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
                sage: D2 = DifferenceOperators(base_ring=P, I=I2, shift=2*h, variable_prefix='w', difference_operator_prefix='D', additional_variables=('t','s'))
                sage: w = D2.functional_variables(); d = D2.difference_operators()
                sage: element1 = w[1,0]*d[0,0]**(-1); element1
                w_1_0*D_0_0^-1
                sage: element2 = w[0,0]*d[2,0]; element2
                w_0_0*D_2_0
                sage: element1*element2
                (w_0_0*w_1_0+(-2*h)*w_1_0)*D_0_0^-1*D_2_0
            """
            parent = self.parent()
            def _shift_degrees(el,shift):
                new_element_dict = {}
                for monom,coeff in el.monomial_coefficients().items():
                    new_element_dict[monom*shift] = coeff
                return parent._from_dict(new_element_dict)
            
            new_element = parent.zero()
            for monom_other,coeff_other in other.monomial_coefficients().items():
                new_element += _shift_degrees(self._lmul_(coeff_other),monom_other)
            return new_element
        
def _to_str(a):
    r"""
    Auxiliary function: return string ``a[1]_a[2]_..._a[k]`` where a=(a[1],...,a[k]) (mathematical notation, not necessarily tuple).

    EXAMPLES::

        sage: _to_str('a')
        'a'
        sage: _to_str(['a','b','c'])
        'a_b_c'
        sage: _to_str([1,2,3])
        '1_2_3'
    """
    a = MonadTuple(a)
    return '_'.join(str(i) for i in a)

    
def MonadTuple(a):
    r"""
    Auxiliary function: return
        - ``tuple(a)`` if ``a`` is iterable
        - ``(a,)`` : `tuple` otherwise.

    EXAMPLES::

        sage: MonadTuple((1,2,3))
        (1, 2, 3)
        sage: MonadTuple([1,2,3])
        (1, 2, 3)
        sage: MonadTuple(1)
        (1,)
        sage: MonadTuple(1) == MonadTuple([1])
        True
    """
    try:
        return tuple(a)
    except:
        return (a,)