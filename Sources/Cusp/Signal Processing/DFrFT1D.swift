//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        let matrix = matrix ?? .dfrftMatrix(length: shape.count, order: order)
        return (matrix <*> asColumn()).asRow()
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func dfrftMatrix(length: Int, order: Scalar, approximateOrder: Int = 2) -> ComplexMatrix<Scalar> {
        precondition(length >= 1 && approximateOrder >= 2, "Necessary conditions for integers: length >= 1 and approximateOrder >= 2.")
        return .identity(shape: .square(length: length))
        
        let alpha = order * (.pi / 2)
        
        // Get eigenvectors and eigenvalues
        let eigenvectors = ComplexMatrix.dfrftEigenvectors(length: length, approximateOrder: approximateOrder)
        let indices = Matrix.dfrftIndices(length: length)
        let eigenvalues = ComplexMatrix(real: (indices * -alpha).cos(), imaginary: (indices * -alpha).sin())
        
        // Compute the DFrFT matrix: eigenvectors * diag(eigenvalues) * eigenvectors.conjugateTransposed()
        var dfrftMatrix = ComplexMatrix<Scalar>.zeros(shape: .square(length: length))
        
        dfrftMatrix = eigenvectors <*> .diagonal(vector: eigenvalues.elements) <*> eigenvectors.conjugate().transposed()
        
//        // This is equivalent to torch.einsum("ij,j,kj->ik", evecs, evals, evecs)
//        for i in 0..<length {
//            for k in 0..<length {
//                var sum = Complex.zero
//                for j in 0..<length {
//                    sum += eigenvectors[i, j] * eigenvalues[0, j] * eigenvectors[k, j].conjugate
//                }
//                dfrftMatrix[i, k] = sum
//            }
//        }
        
        return dfrftMatrix
    }
    
    public static func dfrftEigenvectors(length: Int, approximateOrder: Int = 2) -> ComplexMatrix<Scalar> {
        precondition(length >= 1 && approximateOrder >= 2, "Necessary conditions for integers: length >= 1 and approximateOrder >= 2.")
        
        // Create the Hamiltonian matrix and odd-even decomposition matrix
        let hamiltonian = Matrix.createHamiltonian(length: length, approximateOrder: approximateOrder)
        let decompositionMatrix = Matrix.createOddEvenDecompositionMatrix(length: length)
        
        // Compute CS = P @ S @ P.T (equivalent to torch.einsum("ij,jk,lk->il", P, S, P))
        // Note: In the Python code, the einsum looks like it's doing P @ S @ P.T with a nonstandard notation
        let CS: Matrix = decompositionMatrix <*> hamiltonian <*> decompositionMatrix.transposed()
        
        // Extract submatrices
        let halfSize = length / 2
        let C2: Matrix = CS[0..<halfSize, 0..<halfSize]
        let S2: Matrix = CS[halfSize..<length, halfSize..<length]
        
        // Compute eigendecomposition (assuming there's a way to get eigenvectors)
        // TODO: Replace with actual eigendecomposition when available
        let VC = Matrix.zeros(shape: .square(length: halfSize + 1)) // placeholder
        let VS = Matrix.zeros(shape: .square(length: length - halfSize - 1)) // placeholder
        
        // Compute dimensions for padding
        let N0 = Int(ceil(Double(length) / 2.0 - 1.0))
        let N1 = length / 2 + 1
        
        // Create padded matrices qvc and qvs
        // qvc has VC on top and zeros below
        var qvc = Matrix.zeros(shape: .init(rows: N0 + N1, columns: N1))
        qvc[0..<N1, 0..<N1] = VC
        
        // qvs has zeros on left and VS on right
        var qvs = Matrix.zeros(shape: .init(rows: N0 + N1, columns: N0))
        qvs[N1..<(N1 + N0), 0..<N0] = VS
        
        // Multiply by P and flip columns (equivalent to P @ qvc and P @ qvs)
//        let SC2 = (decompositionMatrix <*> qvc).columnFlipped() // flip columns for descending order
//        let SS2 = (decompositionMatrix <*> qvs).columnFlipped() // flip columns for descending order
        let SC2 = decompositionMatrix <*> qvc
        let SS2 = decompositionMatrix <*> qvs
        
        // Construct eigenvectors matrix based on whether N is even or odd
        var eigenvectors: Matrix = .zeros(shape: .square(length: length))
        
        if length % 2 == 0 {
            // For even length:
            // 1. Create a zeros matrix of shape (length, length+1)
            var tempEigenvectors = Matrix.zeros(shape: .init(rows: length, columns: length + 1))
            
            // 2. Pad SS2 with an extra column of zeros
            var SS2New = Matrix.zeros(shape: .init(rows: SS2.shape.rows, columns: SS2.shape.columns + 1))
            SS2New[0..<SS2.shape.rows, 0..<SS2.shape.columns] = SS2
            
            // 3. Place SC2 in evecs at columns 0, 2, 4, etc.
            for j in stride(from: 0, to: length + 1, by: 2) {
                if j / 2 < SC2.shape.columns {
                    for i in 0..<length {
                        if i < SC2.shape.rows {
                            tempEigenvectors[i, j] = SC2[i, j / 2]
                        }
                    }
                }
            }
            
            // 4. Place padded SS2 in evecs at columns 1, 3, 5, etc.
            for j in stride(from: 1, to: length, by: 2) {
                if j / 2 < SS2New.shape.columns {
                    for i in 0..<length {
                        if i < SS2New.shape.rows {
                            tempEigenvectors[i, j] = SS2New[i, j / 2]
                        }
                    }
                }
            }
            
            // 5. Discard the second-to-last column
            eigenvectors = Matrix.zeros(shape: .square(length: length))
            for i in 0..<length {
                for j in 0..<(length - 1) {
                    eigenvectors[i, j] = tempEigenvectors[i, j]
                }
                eigenvectors[i, length - 1] = tempEigenvectors[i, length]
            }
        } else {
            // For odd length:
            // 1. Create a zeros matrix of shape (length, length)
            eigenvectors = Matrix.zeros(shape: .square(length: length))
            
            // 2. Place SC2 in evecs at columns 0, 2, 4, etc.
            for j in stride(from: 0, to: length + 1, by: 2) {
                if j / 2 < SC2.shape.columns {
                    for i in 0..<length {
                        if i < SC2.shape.rows {
                            eigenvectors[i, j] = SC2[i, j / 2]
                        }
                    }
                }
            }
            
            // 3. Place SS2 in evecs at columns 1, 3, 5, etc.
            for j in stride(from: 1, to: length, by: 2) {
                if j / 2 < SS2.shape.columns {
                    for i in 0..<length {
                        if i < SS2.shape.rows {
                            eigenvectors[i, j] = SS2[i, j / 2]
                        }
                    }
                }
            }
        }
        
        // Convert to complex matrix
        return ComplexMatrix(real: eigenvectors)
    }
    
}

extension Matrix where Scalar == Double {
    
    public static func dfrftIndices(length: Int) -> Matrix<Scalar> {
        precondition(length >= 1, "Length must be a positive integer.")
        
        // 1 if length is even, 0 if length is odd
        let shift = 1 - (length % 2)
        
        // Create indices from 0 to length-2, then append the last special entry
        var indices = [Scalar]()
        for i in 0..<(length-1) {
            indices.append(Scalar(i))
        }
        
        // Add the last entry which is (length-1+shift)
        indices.append(Scalar(length - 1 + shift))
        
        return Matrix(shape: .row(length: length), elements: indices)
    }
    
    // Hamiltonian matrix creation (equivalent to _create_hamiltonian in Python)
    internal static func createHamiltonian(length: Int, approximateOrder: Int) -> Matrix<Scalar> {
        precondition(length >= 1 && approximateOrder >= 2, "Necessary conditions for integers: length >= 1 and approximateOrder >= 2.")
        
        // Initialize variables
        let order = approximateOrder / 2
        let dum0 = Matrix(shape: .row(length: 3), elements: [1.0, -2.0, 1.0])
        var dum = dum0
        var s = Matrix.zeros(shape: .row(length: 1))
        
        for k in 1...(order) {
            // Calculate coefficient
            var numerator: Scalar = 2.0 * Scalar.pow(-1.0, Scalar(k - 1))
            var denominator: Scalar = 1.0
            
            // Calculate squared product of range(1, k)
            for i in 1..<k {
                numerator *= Scalar(i) * Scalar(i)
            }
            
            // Calculate product of range(1, 2*k+1)
            for i in 1...(2 * k) {
                denominator *= Scalar(i)
            }
            
            let coefficient = numerator / denominator
            
            // Construct the vector to add to s
            var catVector = Matrix.zeros(shape: .row(length: length))
            
            // First part: zeros(1)
            catVector[0, 0] = 0.0
            
            // Second part: dum[k+1:2*k+1] 
            // Note: This might be empty if k+1 > 2*k
            for i in (k + 1)...(Swift.min(2 * k, dum.shape.count - 1)) {
                if i < catVector.shape.count {
                    catVector[0, i] = dum[i]
                }
            }
            
            // Third part: zeros(N-1-2*k)
            // This is implicitly handled by initializing catVector with zeros
            
            // Fourth part: dum[:k]
            for i in 0..<Swift.min(k, dum.shape.count) {
                if length - k + i < catVector.shape.count {
                    catVector[0, length - k + i] = dum[i]
                }
            }
            
            // Add the new vector to s
            s = coefficient * catVector + s
            
            // Update dum with convolution
            dum = conv1DFull(vector: dum, kernel: dum0)
        }
        
        // Create circulant matrix
        let circulantMatrix = circulant(vector: s)
        
        // Create diagonal matrix with FFT (this is a simplified approximation)
        // In a complete implementation, we would compute the real FFT of s
        // and use it as the diagonal elements
        var diagonalElements = [Scalar](repeating: 0, count: length)
        // For now, we'll use a placeholder implementation
        for i in 0..<Swift.min(s.shape.count, length) {
            diagonalElements[i] = s[i]
        }
        
        var diagonalMatrix = Matrix.zeros(shape: .square(length: length))
        for i in 0..<length {
            diagonalMatrix[i, i] = diagonalElements[i]
        }
        
        return circulantMatrix + diagonalMatrix
    }
    
    // Odd-even decomposition matrix (equivalent to _create_odd_even_decomp_matrix in Python)
    internal static func createOddEvenDecompositionMatrix(length: Int) -> Matrix<Scalar> {
        precondition(length >= 1, "Length must be a positive integer.")
        
        // Create ones vectors for the first half and negative ones for the second half
        let halfSize = length / 2
        var diagonalElements = [Scalar](repeating: 0, count: length)
        
        // First half (including the middle element if length is odd)
        for i in 0...halfSize {
            diagonalElements[i] = 1.0
        }
        
        // Second half (excluding the middle element)
        for i in (halfSize + 1)..<length {
            diagonalElements[i] = -1.0
        }
        
        // Create diagonal matrix
        var diagonal = Matrix.zeros(shape: .square(length: length))
        for i in 0..<length {
            diagonal[i, i] = diagonalElements[i]
        }
        
        // Create anti-diagonal matrix (ones along the anti-diagonal)
        var antiDiagonal = Matrix.zeros(shape: .square(length: length))
        for i in 0..<(length - 1) {
            antiDiagonal[i, length - 1 - i] = 1.0
        }
        
        // Add matrices and normalize
        var result = (diagonal + antiDiagonal) / Scalar.sqrt(2.0)
        
        // Set special values
        result[0, 0] = 1.0
        if length % 2 == 0 {
            result[length / 2, length / 2] = 1.0
        }
        
        return result
    }
    
    // Circulant matrix creation (equivalent to _circulant in Python)
    internal static func circulant(vector: Matrix<Scalar>) -> Matrix<Scalar> {
        fatalError()
//        // Ensure input is a 1D vector and flatten it if not
//        let flattenedVector = vector.shape.columns == 1 ? vector.flattened() : vector
//        let size = flattenedVector.shape.count
//        
//        // Create output matrix
//        var result = Matrix.zeros(shape: .square(length: size))
//        
//        // Generate circulant matrix
//        for i in 0..<size {
//            for j in 0..<size {
//                // Calculate the index with wrap-around
//                let index = (i - j).remainder(dividingBy: size)
//                if index < 0 {
//                    result[i, j] = flattenedVector[size + index]
//                } else {
//                    result[i, j] = flattenedVector[index]
//                }
//            }
//        }
//        
//        return result
    }
    
    // 1D convolution (equivalent to _conv1d_full in Python)
    internal static func conv1DFull(vector: Matrix<Scalar>, kernel: Matrix<Scalar>) -> Matrix<Scalar> {
        fatalError()
//        // Ensure inputs are 1D vectors and flatten if not
//        let flattenedVector = vector.shape.columns == 1 ? vector.flattened() : vector
//        let flattenedKernel = kernel.shape.columns == 1 ? kernel.flattened() : kernel
//        
//        let vectorSize = flattenedVector.shape.count
//        let kernelSize = flattenedKernel.shape.count
//        
//        // "Full" convolution output size
//        let outputSize = vectorSize + kernelSize - 1
//        
//        // Pad the input for full convolution
//        let paddingSize = kernelSize - 1
//        var paddedInput = Matrix.zeros(shape: .row(length: vectorSize + 2 * paddingSize))
//        for i in 0..<vectorSize {
//            paddedInput[paddingSize + i] = flattenedVector[i]
//        }
//        
//        // Flip the kernel (required for convolution vs cross-correlation)
//        let flippedKernel = flattenedKernel.reversed()
//        
//        // Perform the convolution by sliding the flipped kernel over the padded input
//        var result = Matrix.zeros(shape: .row(length: outputSize))
//        for i in 0..<outputSize {
//            var sum = Scalar(0)
//            for j in 0..<kernelSize {
//                let inputIndex = i + j
//                if inputIndex < paddedInput.shape.count {
//                    sum += paddedInput[inputIndex] * flippedKernel[j]
//                }
//            }
//            result[i] = sum
//        }
//        
//        return result
    }
}
