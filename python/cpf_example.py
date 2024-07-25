from mHT.CPF import cpf_fast, cpf_accurate

#example parameters for the cpf functions
x=1; #dimensionless
y=1; #dimensionless

print("the output of the cpf_accurate function")
print(cpf_accurate(x,y))

print("the output of the cpf_fast function")
print(cpf_fast(x,y))
