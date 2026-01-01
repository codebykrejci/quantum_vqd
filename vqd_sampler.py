from scipy.optimize import minimize

import numpy as np

import datetime

import h5py

import time



from qiskit.circuit import QuantumCircuit, ParameterVector


import config_bg, config_CuO2, config_Si


import cmp




config = config_CuO2     # for bilayer graphene or silicon just change config to config_bg or config_Si

Nshots = config.NSHOTS

simulator = cmp.SIMULATOR





def calculate_overlaps(ansatz: QuantumCircuit, 

                       previous_circuits: list[tuple[QuantumCircuit, np.ndarray]], 

                       variational_parameters: dict) -> np.ndarray:

    """

    Calculates the overlap terms between the state vector of a parameterized quantum circuit

    (ansatz) and a state vectors from previously optimized quantum circuits using the inverse circuit method.



    :param ansatz: The quantum circuit representing the current ansatz.

    :type ansatz: qiskit.circuit.QuantumCircuit

    :param previous_circuits_and_params: A list of tuples, where each tuple contains a

                                         QuantumCircuit and its optimized parameters (NumPy array)

                                         for the previous states.

    :type previous_circuits_and_params: list[tuple[qiskit.circuit.QuantumCircuit, numpy.ndarray]]

    :param variational_parameters: A dictionary mapping parameters of the ansatz

                                   to their numerical values for the current state.

    :type variational_parameters: dict

    :returns: A NumPy array containing the squared absolute overlaps

              (|<current_state|previous_state>|^2) between the current ansatz

              state and each of the previous circuit states.

    :rtype: numpy.ndarray

    """

    overlaps = []

    current_state = ansatz.assign_parameters(variational_parameters)



    for prev_circuit in previous_circuits:

        curr = current_state.copy()

        prev = prev_circuit.copy()



        fidelity_circuit = curr.compose(prev.inverse())

        fidelity_circuit.measure_all()



        result = simulator.run(fidelity_circuit, shots=Nshots).result()

        zero_bitstring = "0" * fidelity_circuit.num_qubits

        overlap = result.get_counts().get(zero_bitstring, 0) / Nshots



        overlaps.append(overlap)



    return np.array(overlaps)



def cost_func_vqd(parameters: np.ndarray, 

                  ansatz: QuantumCircuit,

                  num_qubits: int,

                  hamiltonian: np.ndarray, 

                  previous_states: list[tuple[QuantumCircuit, np.ndarray]], 

                  betas: np.ndarray) -> float:

    """

    Cost function for Variational Quantum Deflation (VQD) using measurement statistics



    This function calculates the energy expectation value of the qubit Hamiltonian with respect

    to the given ansatz, and then adds a penalty term based on overlaps with previously

    found eigenstates. 



    :param parameters: A 1D NumPy array containing the current numerical values for the variational

                       parameters of the ansatz circuit. These values are iteratively adjusted

                       by the optimization algorithm to minimize the cost function.

    :type parameters: numpy.ndarray

    :param ansatz: The quantum circuit representing the current ansatz.

                   The circuit's parameters are optimized during the VQD process.

    :type ansatz: qiskit.circuit.QuantumCircuit

    :param num_qubits: Number of qubits.

    :type num_qubits: integer

    :param hamiltonian: The Hamiltonian operator for which the energy expectation value

                        is to be calculated. 

    :type hamiltonian: numpy.ndarray

    :param previous_states: A list of previously optimized quantum states. Each element in this list

                        is a tuple containing a `QuantumCircuit` (with the same ansatz structure)

                        and a NumPy array of its optimal numerical parameters.

    :type prev_states: list[tuple[qiskit.circuit.QuantumCircuit, numpy.ndarray]]

    :param betas: A NumPy array containing the deflation coefficients.

                  These coefficients determine the magnitude of the penalty for the overlap terms.

                  They are crucial for ensuring the VQD algorithm finds distinct eigenstates.

    :type betas: numpy.ndarray

    :returns: The total cost for the VQD optimization. This includes the energy

              expectation value and the deflation penalty, which the optimizer aims to minimize.

    :rtype: float

    """

    bound_circuit = ansatz.assign_parameters(parameters)

    probabilities, amplitudes_abs, nonzero_keys, nonzero_pairs = cmp.amplitudes(bound_circuit, Nshots)

    N = bound_circuit.num_qubits

    cost = 0.0

    

    for j in range(N):

        epsilon = hamiltonian[j, j].real

        if abs(epsilon) > 1e-5:

            cost += epsilon * probabilities[j]

    

    

    if len(nonzero_keys) > 1:

        exp_xx = cmp.expectation_values_xx(bound_circuit, Nshots, nonzero_keys)

        exp_xy = cmp.expectation_values_xy(bound_circuit, Nshots, nonzero_keys)

        complex_corr = cmp.C(exp_xx, exp_xy)

        

        for (j, l) in nonzero_pairs:

            if (j, l) in complex_corr:

                cost += (hamiltonian[j, l] * complex_corr[(j, l)]).real

            else:

                for k_mid in range(j + 1, l):

                    if (j, k_mid) in complex_corr and (k_mid, l) in complex_corr:

                        cost += (hamiltonian[j, l] * complex_corr[(j, k_mid)] *complex_corr[(k_mid, l)] /(2 * probabilities[k_mid])).real

                        

        

    if previous_states:

        overlaps = calculate_overlaps(ansatz, [state.assign_parameters(params) for state, params in previous_states], parameters)

        cost += np.sum([np.real(betas[i] * overlaps[i]) for i in range(len(overlaps))])

    

    return cost



class QuantumBandStructureCalculator:

    def __init__(self):

        

        self.num_qubits = config.NUM_QUBITS

        self.neighbors = config.NEIGHBORS

        self.hpar = config.PARAMS

        self.calc_method = config.calc_method

        self.max_iter = config.max_iter

        self.bootstrapping = config.bootstrapping

        self.num_states = config.num_states

        self.path_q = config.path_q

        self.exact_path_q = config.path_exact



        self.circuit = self._build_ansatz_circuit()

    

    def _build_ansatz_circuit(self) -> QuantumCircuit:

        """

        Builds the variational ansatz circuit based on the number of qubits.



        This method constructs a quantum circuit with a specific parameterized structure. 

        It initializes a ParameterVector to manage the circuit's variational parameters.



        :returns: The constructed Qiskit QuantumCircuit representing the ansatz.

        :rtype: qiskit.circuit.QuantumCircuit

        """

            

        circuit = QuantumCircuit(self.num_qubits)

        params = ParameterVector('theta', 2*(self.num_qubits-1))  



        circuit.x(0)



        # Helper function for applying A-gate

        def _apply_ops(circuit, qubit1, qubit2, param1, param2):

            circuit.cx(qubit1, qubit2)

            circuit.rz(-param1 - np.pi, qubit1)

            circuit.ry(-param2 - np.pi/2, qubit1)

            circuit.cx(qubit2, qubit1)

            circuit.ry(param2 + np.pi/2, qubit1)

            circuit.rz(param1 + np.pi, qubit1)

            circuit.cx(qubit1, qubit2)



        for i in range(self.num_qubits - 1):  

            _apply_ops(circuit, i, i+1, params[2*i], params[2*i+1])



        return circuit

    

    def energies_q(self,k):

        """

        Computes the eigenvalues (energies) of a given Hamiltonian

        at point in k-space using the Variational Quantum Deflation (VQD) method.



        :param k: A k-vector (wave vector) used in building the Hamiltonian.

        :type k: numpy.ndarray



        :returns: A tuple containing:

            - eigenvalues_sorted (numpy.ndarray): Sorted list of approximated eigenvalues (energies).

            - n_fun_sorted (numpy.ndarray): Sorted list of the number of function evaluations per eigenvalue.

            - opt_parameters_sorted (list[numpy.ndarray]): List of optimal parameters (one array per state), sorted by eigenvalue.

            - time_duration_sorted (list[float]): List of durations for each minimization step, sorted by eigenvalue.

        :rtype: tuple[numpy.ndarray, numpy.ndarray, list[numpy.ndarray], list[float]]



        :notes:

            - Bootstrapping is supported: if `self.bootstrapping` is enabled and `self.step > 0`, 

            previously calculated optimal parameters are reused as starting points.

        """



        x0 = np.asarray([np.random.uniform(0, 2 * np.pi, 2*(self.num_qubits-1))] * self.num_states) 

        prev_states = [] 

        eigenvalues = []

        opt_parameters = []

        n_fun = [] 

        time_duration = []

        class_ham, class_eignvls = config.hamiltonian(k, *self.hpar)

        upper_bound = np.sum(np.abs(np.triu(class_ham, k=0)))

        betas = np.asarray([upper_bound * 10] * (self.num_states - 1))





        for n_step in range(self.num_states):

            print(f'Step {n_step + 1}: Starting minimization')

            start_time_minimize = time.time()



            initial_parameters = x0[n_step]

            if (self.step > 0) & (self.bootstrapping):

                initial_parameters = self.calc_optimal_param[n_step] 



            result = minimize(cost_func_vqd,

                            x0=initial_parameters,

                            args=(self.circuit, self.num_qubits, class_ham, prev_states, betas),

                            method=self.calc_method,

                            options={'maxiter': self.max_iter, 'initial_tr_radius': 0.5, 'final_tr_radius': 0.01},

                            tol=1e-3)

            

            end_time_minimize = time.time()

            duration_minimize = end_time_minimize - start_time_minimize

            print(f'Minimization done, duration: {duration_minimize:.4f} s')

            time_duration.append(duration_minimize)



            opt_parameters.append(result.x)

            eigenvalues.append(result.fun)

            n_fun.append(result.nfev)

            prev_states.append((self.circuit, result.x)) 

            



        sorted_ind = np.argsort(eigenvalues)

        eigenvalues_sorted = np.array(eigenvalues)[sorted_ind] 

        n_fun_sorted = np.array(n_fun)[sorted_ind] 

        opt_parameters_sorted = [opt_parameters[i] for i in sorted_ind]

        time_duration_sorted = [time_duration[i] for i in sorted_ind]

        print("Energies:", eigenvalues_sorted)



        return eigenvalues_sorted, n_fun_sorted, opt_parameters_sorted, time_duration_sorted

    

    def run_calculation(self):

        """

        Runs the full bandstructure calculation over all k-points in the defined path,

        stores results to an HDF5 file, and computes reference (exact) energies for comparison using the exact diagonalization method.



        :returns: None



        :effects:

            - Writes results into a file named `results_YYYYMMDD-HHMMSS.h5`.

            - Updates internal state (`self.step`, `self.calc_optimal_param`) for bootstrapping.

            - Prints time durations and progress to stdout.



        :file structure:

            - Scalar attributes:

                - `calc_method`, `max_iter`, `bootstrapping`, `num_states`, `lattice_constant`, `path_print`

            - For each k-point:

                - Group `k-point index {i}` with datasets:

                    - `eigenvalues` – sorted list of approximated eigenenergies

                    - `n_fun` – number of function evaluations for each state

                    - `optimal_params` – optimal parameters for each state

                    - `minimize_time` – duration of each optimization step  

            - `exact_values` group:

                - `eigenvalues` – list of classically computed eigenvalues

                - `path` – the corresponding path used for exact values



        """



        start_time_script = time.time()

        self.step = 0

        timename = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')



        with h5py.File(f'results_{timename}.h5', 'w') as f:

            f.create_dataset('calc_method', data = self.calc_method, dtype=h5py.string_dtype(encoding='utf-8'))

            f.create_dataset('max_iter', data = self.max_iter)

            f.create_dataset('bootstrapping', data = self.bootstrapping)

            f.create_dataset('num_states', data = self.num_states)

            f.create_dataset('lattice_constant', data = config.LATTICE_CONSTANT)

            f.create_dataset('path_print', data = config.path_q_plot)

            f.create_dataset('labels_position', data = config.labels_position)

            f.create_dataset('labels_name', data = config.labels_name, dtype=h5py.string_dtype(encoding='utf-8'))



            print('Calculating exact values...')

            exact_values = []

            for k_ in self.exact_path_q:

                class_ham, class_eignvls = config.hamiltonian(k_, *self.hpar)

                exact_values.append(class_eignvls)



            exact_group = f.create_group('exact_values')

            exact_group.create_dataset('eigenvalues', data = exact_values)

            exact_group.create_dataset('path', data = config.path_exact_plot)

            print('Done.')

            print('Starting with quantum computation.')

            

            start_calc_time = time.time()



            for k_ in range(len(self.path_q)):

                k_group = f.create_group(f'k-point index {k_}')

                print(f'Calculating {k_+1}/{len(self.path_q)}')

                result = self.energies_q(self.path_q[k_])



                k_group.create_dataset('eigenvalues', data=result[0])

                k_group.create_dataset('n_fun', data=result[1])

                k_group.create_dataset('optimal_params', data=result[2])

                k_group.create_dataset('minimize_time', data=result[3]) 



                f.flush() 



                self.calc_optimal_param = result[2]

                self.step += 1



            end_calc_time_script = time.time()

            duration_calc = end_calc_time_script - start_calc_time

            print(f'Total time duration for calculated values: {duration_calc:.4f} s')

        

            f.create_dataset('calculated_values_duration', data=duration_calc)





        end_time_script = time.time()

        duration_calc = end_time_script - start_time_script

        print(f'Total time duration: {duration_calc:.4f} s')





if __name__ == '__main__':

    qbsc = QuantumBandStructureCalculator()

    qbsc.run_calculation()

    

