from mHT.CPF import cpf_fast, cpf_accurate

# Example parameters for the cpf functions
x=1 # Dimensionless
y=1 # Dimensionless

print("The output of the cpf_accurate function:")
print(cpf_accurate(x,y))

print("The output of the cpf_fast function:")
print(cpf_fast(x,y))
