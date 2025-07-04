def extract_elastic_constants_from_file(file_path):
    pattern = "Elastic Stiffness Constants Cij (GPa)"
    with open(file_path, 'r') as file:
        lines = file.readlines()
        start_line = None

        for i, line in enumerate(lines):
            line = line.strip()
            if pattern in line:
                start_line = i + 3
                break
        if start_line is not None:
            extracted_data = lines[start_line:start_line + 6]
            data_list = []
            for line in extracted_data:
                data_list.extend(float(number) for number in line.split())
            
            if len(data_list) == 36:
                c11, c12, c13 = data_list[0], data_list[1], data_list[2]
                c14, c15, c16 = data_list[3], data_list[4], data_list[5]
                c21, c22, c23 = data_list[6], data_list[7], data_list[8]
                c24, c25, c26 = data_list[9], data_list[10], data_list[11]
                c31, c32, c33 = data_list[12], data_list[13], data_list[14]
                c34, c35, c36 = data_list[15], data_list[16], data_list[17]
                c41, c42, c43 = data_list[18], data_list[19], data_list[20]
                c44, c45, c46 = data_list[21], data_list[22], data_list[23]
                c51, c52, c53 = data_list[24], data_list[25], data_list[26]
                c54, c55, c56 = data_list[27], data_list[28], data_list[29]
                c61, c62, c63 = data_list[30], data_list[31], data_list[32]
                c64, c65, c66 = data_list[33], data_list[34], data_list[35]
                return (
                    c11, c12, c13, c14, c15, c16,
                    c21, c22, c23, c24, c25, c26,
                    c31, c32, c33, c34, c35, c36,
                    c41, c42, c43, c44, c45, c46,
                    c51, c52, c53, c54, c55, c56,
                    c61, c62, c63, c64, c65, c66
                )

            else:
                return "Data list and constants list lengths do not match."
        else:
            return "'Elastic Stiffness Constants Cij (GPa)' not found in the file."

def extract_other_constants_from_file(file_path):
    pattern = "Elastic constants for polycrystalline material (GPa)"
    with open(file_path, 'r') as file:
        lines = file.readlines()
        start_line = None

        for i, line in enumerate(lines):
                line = line.strip()
                if pattern in line:
                    start_line = i + 4
                    break
        if start_line is not None:
            extracted_data = lines[start_line:start_line + 5]
            data_list = []
            for line in extracted_data:
                values = line.split(':', 1)
                data_list.extend(float(number) for number in values[1].split())
            # bulk_modulus_voigt, bulk_modulus_reuss, bulk_modulus = data_list[0], data_list[1], data_list[2]
            # shear_modulus_voigt, shear_modulus_reuss, shear_modulus = data_list[3], data_list[4], data_list[5]
            # lame_lambda_voigt, lame_lambda_reuss, lame_lambda = data_list[6], data_list[7], data_list[8]
            # young_modulus_voigt, young_modulus_resuss, young_modulus = data[9], data[10], data[11]
            # poisson_ratio_voigt, possion_ratio_resuss, poisson_ratio = data[12], data[13], data[14]
            constant_dictonary = {'bulk_modulus_voigt': data_list[0], 
                                   'bulk_modulus_reuss': data_list[1],
                                   'bulk_modulus': data_list[2],
                                   'shear_modulus_voigt': data_list[3], 
                                   'shear_modulus_reuss': data_list[4], 
                                   'shear_modulus': data_list[5],
                                   'lame_lambda_voigt': data_list[6], 
                                   'lame_lambda_reuss': data_list[7],
                                   'lame_lambda': data_list[8],
                                   'young_modulus_voigt': data_list[9], 
                                   'young_modulus_reuss': data_list[10],
                                   'young_modulus' : data_list[11], 
                                   'poisson_ratio_voigt': data_list[12], 
                                   'poisson_ratio_reuss': data_list[13],
                                   'poisson_ratio': data_list[14]}
            return constant_dictonary
        else:
            return "'Elastic constants for polycrystalline material (GPa)' not found in the file."

def extract_cell_volume_from_file(file_path):
    pattern = "Current cell volume ="
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            cell_volume = float(last_match.split('=')[1].strip().split(' ')[0])
            return cell_volume * 1e-30
        else:
            print("No match found")

def extract_density_from_file(file_path):
    pattern = 'density ='
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            density = float(last_match.split('=')[1].strip().split(' ')[0])
            return density * 1660.54
        else:
            print("No match found")
def extract_compressibility_from_file(file_path):
    pattern = 'Compressibility'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            compressibility = float(last_match.split('=')[1].strip().split(' ')[0])
            return compressibility
        else:
            print("No match found")

def extract_elastic_debye_temperature_from_file(file_path):
    pattern = 'Elastic Debye temperature'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            elastic_debye_temperature = float(last_match.split('=')[1].strip().split(' ')[0])
            return elastic_debye_temperature
        else:
            print("No match found")

def extract_averaged_sound_velocity_from_file(file_path):
    pattern = 'Averaged sound velocity'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            averaged_sound_velocity = float(last_match.split('=')[1].strip().split(' ')[0])
            return averaged_sound_velocity
        else:
            print("No match found")

def extract_number_of_ions_from_file(file_path):
    pattern = 'number of ions'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            number_of_ions = float(last_match.split('=')[1].strip())
            return number_of_ions 
        else:
            print("No match found")

def extract_number_of_species_from_file(file_path):
    pattern = 'number of species'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            number_of_species = float(last_match.split('=')[1].strip())
            return number_of_species 
        else:
            print("No match found")
            
def extract_max_number_of_species_formula_unit_from_file(file_path):
    pattern = 'Max number of any one species'
    last_match = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if pattern in line:
                last_match = line.strip()
        if last_match:
            formula_factor = float(last_match.split('=')[1].strip())
            return formula_factor 
        else:
            print("No match found")
                
def print_non_zero_value_of_Cij(elastic_constants):
    count = 0
    for i in range(6):
        for j in range(6):
            if elastic_constants[count] != 0:
                print(f"c{i+1}{j+1} {' ' * (20 - len(f'c{i+1}{j+1}'))}: {elastic_constants[count]}")
            count+=1

if __name__ == "__main__":
    file_path = r'F:\elastic constant calculations\example.txt'
    file_path1 = 'MgSiP2.castep'
    

    c11, c12, c13, c14, c15, c16 = extract_elastic_constants_from_file(file_path)[:6]
    c21, c22, c23, c24, c25, c26 = extract_elastic_constants_from_file(file_path)[6:12]
    c31, c32, c33, c34, c35, c36 = extract_elastic_constants_from_file(file_path)[12:18]
    c41, c42, c43, c44, c45, c46 = extract_elastic_constants_from_file(file_path)[18:24]
    c51, c52, c53, c54, c55, c56 = extract_elastic_constants_from_file(file_path)[24:30]
    c61, c62, c63, c64, c65, c66 = extract_elastic_constants_from_file(file_path)[30:36]


    other_constants = extract_other_constants_from_file(file_path)
    bulk_modulus_voigt, bulk_modulus_reuss = other_constants['bulk_modulus_voigt'], other_constants['bulk_modulus_reuss']
    bulk_modulus = other_constants['bulk_modulus']
    shear_modulus_voigt, shear_modulus_reuss = other_constants['shear_modulus_voigt'], other_constants['shear_modulus_reuss']
    shear_modulus = other_constants['shear_modulus']
    lame_lambda_voigt, lame_lambda_reuss = other_constants['lame_lambda_voigt'], other_constants['lame_lambda_reuss']
    lame_lambda = other_constants['lame_lambda']
    young_modulus_voigt, young_modulus_resuss = other_constants['young_modulus_voigt'], other_constants['young_modulus_reuss']
    young_modulus = other_constants['young_modulus']
    poisson_ratio_voigt, possion_ratio_resuss = other_constants['poisson_ratio_voigt'], other_constants['poisson_ratio_reuss'] 
    poisson_ratio = other_constants['poisson_ratio']

    # print(poisson_ratio)
    # print_non_zero_value_of_Cij(extract_elastic_constants_from_file(file_path))
    cell_volume = extract_number_of_ions_from_file(file_path1)
    print(cell_volume)
   
    
