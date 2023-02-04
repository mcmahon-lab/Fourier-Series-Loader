''' 
This files includes codes for various functions that are needed 
for simulating the FSL method. 
'''

import numpy as np

def Fourier_state(f,m):
    
    r'''
    This function calculates the (m+1)-qubit
    quantum states |c> encoding the first 2^(m+1) 
    Fourier modes of the function f.
    
    Note that the this state is defined in Eq.(1)
    of our paper.
        
    Args:
        f (list): Target state of a function 
            whose Fourier state we need to prepare.
        m (int): Determines the number of Fourier
            modes that we need to determine.
            
    Returns:
        c_state (list): Quantum state of m+1 qubits
            with 2^{m+1} dominant Fourier coefficients.
            
    Raises ValueError:
        If m is not less than total number of qubits.
    '''
    
    if len(f) < 2**(m+1):
        raise ValueError("m should be less than total number of qubits.")
    
    # Fourier Coefficients
    c = np.fft.ifft(f)

    # Fourier state:
    c_state = [*c[:2**m],0,*c[len(f)+1-2**m:]]
    c_state = c_state/np.linalg.norm(c_state)
    
    return c_state

def Fourier_state_2d(f,m):
    
    r'''
    This function calculates the 2*(m+1)-qubit
    quantum states |c> encoding the first 2^(2m+2) 
    Fourier modes of a function f of two variables.
    
    Note that the this state is defined in Eq.(C3)
    of our paper.
        
    Args:
        f (list): Target state of a 2d function 
            whose Fourier state we need to prepare.
        m (int): Determines the number of Fourier
            modes that we need to determine.
            
    Returns:
        c_state (list): Quantum state of 2*(m+1) qubits
            with 2^{2m+2} dominant Fourier coefficients.
            
    Raises ValueError:
        If m is not less than total number of qubits.
    '''
    
    N = len(f)
    if N < 2**(m+1):
        raise ValueError("m should be less than total number of qubits per dimension.")
    
    # Fourier Coefficients
    c = np.fft.ifft2(f)

    # Fourier state:
    c_state = []
    
    for k in range(2**m):
        c_state.extend([*c[k][:2**m],0,*c[k][N+1-2**m:]])
    c_state.extend(np.zeros(2**(m+1)))
    for k in range(N+1-2**m,N):
        c_state.extend([*c[k][:2**m],0,*c[k][N+1-2**m:]])
            
    c_state = c_state/np.linalg.norm(c_state)
    
    return c_state


def output_reordering(qiskit_output):
    
    '''
    The output from the Qiskit uses a different qubit
    ordering convention from the convention that we used.
    For example, Qiskit labels the state |01000> as |2>
    and |00010> as |16>. On the other hand, our convention
    demands that |2> should be a label of |00010> and |16> 
    be a label of |01000>. 
    
    This function takes an output from the Qiskit and
    returns the output which is consistent with our
    convention.
    
    Args:
        qiskit_output (list): Output from the simulation
            performed by Qiskit. 
    '''
    
    # We first determine how many qubits there are.
    # Note: dim(qiskit_output) = 2^n
    n = int(np.round(np.log(len(qiskit_output))/np.log(2)))
    
    new_output = []
    
    for k in range(len(qiskit_output)):
        
        bin_k = bin(k)[2:]
        bin_k = (n-len(bin_k))*'0' + bin_k
        bin_k_rev = bin_k[::-1]
        k_rev = int(bin_k_rev,2)
        
        new_output.append(qiskit_output[k_rev])
    
    return new_output   