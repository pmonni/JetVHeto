"""Independent power counting for the two radius prescriptions (SymPy)."""
import sympy as s
a,L,b0,b1,d2,d3,q,h,lam=s.symbols('a L b0 b1 Delta2 Delta3 q h lambda')
w=1-2*lam
t=-2*s.pi*a*b1/b0*s.log(w)/w
A=d2*a**2*L/(4*w**2)
B=d3*a**3*L/(8*w**3)
C=-s.pi*b0*d2*a**3*L*((1-3*lam+2*lam**2)*q+h)/w**3
products=A*(1+t)**2+B*(1+t)**3+C*(1+t)
strict=A*(1+2*t)+B+C
extra=A*t**2+B*(3*t+3*t**2+t**3)+C*t
assert s.simplify(products-strict-extra)==0
fixed_L=s.series(extra.subs(lam,2*s.pi*b0*a*L),a,0,6).removeO().expand()
assert all(fixed_L.coeff(a,k)==0 for k in range(5))
leading=s.pi**2*b1*L**2*(3*d3-8*s.pi*b0*d2*(q+h))
assert s.simplify(fixed_L.coeff(a,5)-leading)==0
fixed_lambda=s.expand(extra.subs(L,lam/(2*s.pi*b0*a)))
assert all(fixed_lambda.coeff(a,k)==0 for k in range(3))
print('Difference starts at a^3 at fixed lambda (N4LL).')
print('At fixed L the difference vanishes through a^4; leading term is:')
print('a^5 *',s.factor(leading))
print('The common luminosity, hard terms and Sudakov cannot lower this first power.')
print('Matching expansion through relative a^3 is therefore unchanged.')
