import regex as re

def convert_powers(equation):
    # Pattern to match terms with ^2 or ^(1/2)
    pattern = r'(\([^\)]+\)|\w+)\^(\d+|\(1/2\))'
    #pattern = r'(\((?:[^\(\)]+|(?R))*\)|\w+)\^(\d+|\(1/2\))'
    
    # Function to replace matched patterns
    def replacer(match):
        base = match.group(1)  # the term before ^
        exponent = match.group(2)  # the power (2 or (1/2))
        
        # Use std::pow for ^2 and std::sqrt for ^(1/2)
        if exponent == "2":
            return f"std::pow({base}, 2)"
        elif exponent == "(1/2)":
            return f"std::sqrt({base})"
        else:
            print("exponent:" + f"{exponent}")
            # For other cases, leave it as is for now
            return f"{base}^{exponent}"
            

    # Substitute all matches in the equation
    converted_equation = re.sub(pattern, replacer, equation)
    
    return converted_equation

# Example usage
equation = "(2*z_k*(Px_0^2 - 2*Px_0*x_kp1 + Py_0^2 - 2*Py_0*y_kp1 + x_kp1^2 + y_kp1^2) + z_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Py_0*y_k + 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 - 2*Px_0^2 - 2*Py_0^2) - Pz_0*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*x_kp1^2 + 2*y_kp1^2) + (z_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 + 2*x_k*x_kp1 + 2*y_k*y_kp1 - 2*x_k^2 - 2*y_k^2))/N_p - (z_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*x_kp1^2 + 2*y_kp1^2))/N_p + (2*Pz_0*(x_k^2 - 2*x_k*x_kp1 + x_kp1^2 + y_k^2 - 2*y_k*y_kp1 + y_kp1^2))/N_p)/(4*(x_k^2*y_kp1^2 + x_kp1^2*y_k^2 + Pz_0^2*(x_k^2 - 2*x_k*x_kp1 + x_kp1^2 + y_k^2 - 2*y_k*y_kp1 + y_kp1^2) + Py_0^2*x_k^2 + Py_0^2*x_kp1^2 + Px_0^2*y_k^2 + Px_0^2*y_kp1^2 + z_kp1^2*(Px_0^2 - 2*Px_0*x_k + Py_0^2 - 2*Py_0*y_k + x_k^2 + y_k^2) + z_k^2*(Px_0^2 - 2*Px_0*x_kp1 + Py_0^2 - 2*Py_0*y_kp1 + x_kp1^2 + y_kp1^2) - 2*Py_0^2*x_k*x_kp1 - 2*Px_0*x_k*y_kp1^2 - 2*Px_0*x_kp1*y_k^2 - 2*Py_0*x_k^2*y_kp1 - 2*Py_0*x_kp1^2*y_k - 2*Px_0^2*y_k*y_kp1 + z_k*z_kp1*(2*Px_0*x_k + 2*Px_0*x_kp1 + 2*Py_0*y_k + 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 - 2*Px_0^2 - 2*Py_0^2) + Pz_0*z_kp1*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 + 2*x_k*x_kp1 + 2*y_k*y_kp1 - 2*x_k^2 - 2*y_k^2) - Pz_0*z_k*(2*Px_0*x_k - 2*Px_0*x_kp1 + 2*Py_0*y_k - 2*Py_0*y_kp1 - 2*x_k*x_kp1 - 2*y_k*y_kp1 + 2*x_kp1^2 + 2*y_kp1^2) + 2*Py_0*x_k*x_kp1*y_k + 2*Py_0*x_k*x_kp1*y_kp1 + 2*Px_0*x_k*y_k*y_kp1 + 2*Px_0*x_kp1*y_k*y_kp1 - 2*x_k*x_kp1*y_k*y_kp1 - 2*Px_0*Py_0*x_k*y_k + 2*Px_0*Py_0*x_k*y_kp1 + 2*Px_0*Py_0*x_kp1*y_k - 2*Px_0*Py_0*x_kp1*y_kp1)^(1/2)) + (2*z_k*(Px_0^2 - 2*Px_0*x_km1 + Py_0^2 - 2*Py_0*y_km1 + x_km1^2 + y_km1^2) + z_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Py_0*y_k + 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 - 2*Px_0^2 - 2*Py_0^2) - Pz_0*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 + 2*x_km1^2 + 2*y_km1^2))/(4*(x_k^2*y_km1^2 + x_km1^2*y_k^2 + Pz_0^2*(x_k^2 - 2*x_k*x_km1 + x_km1^2 + y_k^2 - 2*y_k*y_km1 + y_km1^2) + Py_0^2*x_k^2 + Py_0^2*x_km1^2 + Px_0^2*y_k^2 + Px_0^2*y_km1^2 + z_km1^2*(Px_0^2 - 2*Px_0*x_k + Py_0^2 - 2*Py_0*y_k + x_k^2 + y_k^2) + z_k^2*(Px_0^2 - 2*Px_0*x_km1 + Py_0^2 - 2*Py_0*y_km1 + x_km1^2 + y_km1^2) - 2*Py_0^2*x_k*x_km1 - 2*Px_0*x_k*y_km1^2 - 2*Px_0*x_km1*y_k^2 - 2*Py_0*x_k^2*y_km1 - 2*Py_0*x_km1^2*y_k - 2*Px_0^2*y_k*y_km1 + z_k*z_km1*(2*Px_0*x_k + 2*Px_0*x_km1 + 2*Py_0*y_k + 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 - 2*Px_0^2 - 2*Py_0^2) + Pz_0*z_km1*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 + 2*x_k*x_km1 + 2*y_k*y_km1 - 2*x_k^2 - 2*y_k^2) - Pz_0*z_k*(2*Px_0*x_k - 2*Px_0*x_km1 + 2*Py_0*y_k - 2*Py_0*y_km1 - 2*x_k*x_km1 - 2*y_k*y_km1 + 2*x_km1^2 + 2*y_km1^2) + 2*Py_0*x_k*x_km1*y_k + 2*Py_0*x_k*x_km1*y_km1 + 2*Px_0*x_k*y_k*y_km1 + 2*Px_0*x_km1*y_k*y_km1 - 2*x_k*x_km1*y_k*y_km1 - 2*Px_0*Py_0*x_k*y_k + 2*Px_0*Py_0*x_k*y_km1 + 2*Px_0*Py_0*x_km1*y_k - 2*Px_0*Py_0*x_km1*y_km1)^(1/2))" 
converted_equation = convert_powers(equation)
print("Converted equation:", converted_equation)


def revert_powers(converted_equation):
    # Pattern to match std::pow(base, 2) and std::sqrt(base)
    pattern_pow = r'std::pow\(([^,]+), 2\)'
    pattern_sqrt = r'std::sqrt\(([^)]+)\)'

    # Replace std::pow(base, 2) with base^2
    reverted_equation = re.sub(pattern_pow, r'\1^2', converted_equation)
    
    # Replace std::sqrt(base) with base^(1/2)
    reverted_equation = re.sub(pattern_sqrt, r'\1^(1/2)', reverted_equation)
    
    return reverted_equation

# Example usage
reverted_equation = revert_powers(converted_equation)
print(reverted_equation==equation)

