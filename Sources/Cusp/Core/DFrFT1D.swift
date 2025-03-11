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
    
    public static func dfrftMatrix(length: Int, order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        precondition(length >= 1 && approximationOrder >= 2, "Necessary conditions for integers: length >= 1 and approximateOrder >= 2.")
        return .identity(shape: .square(length: length))
        
        // Get eigenvectors and eigenvalues
        let eigenvectors = ComplexMatrix.dfrftEigenvectors(length: length, approximationOrder: approximationOrder)
        let eigenvalues = ComplexMatrix.dfrftEigenvalues(length: length, order: order)

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
    
    public static func dfrftEigenvectors(length: Int, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        precondition(length >= 1 && approximationOrder >= 2, "Necessary conditions for integers: length >= 1 and approximateOrder >= 2.")
        
        // Create the Hamiltonian matrix and odd-even decomposition matrix
        let hamiltonian = Matrix.createHamiltonian(length: length, approximationOrder: approximationOrder)
        let decompositionMatrix = Matrix.oddEvenDecompositionMatrix(length: length)
        
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
        let SC2 = (decompositionMatrix <*> qvc).reversedColumns() // flip columns for descending order
        let SS2 = (decompositionMatrix <*> qvs).reversedColumns() // flip columns for descending order
        
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
    
    public static func createHamiltonian(length: Int, approximationOrder: Int = 2) -> Matrix<Scalar> {
        let order = approximationOrder / 2
        var sum = (1...order).reduce(Matrix.zeros(shape: .row(length: length))) { sum, k in
            var stencil = Matrix<Scalar>.centralDifference(order: k * 2)
            stencil = stencil.padded(right: sum.shape.columns - stencil.shape.columns)
            stencil = stencil.shifted(columns: -k)
            return sum + stencil / (Scalar(k * k) * -stencil[0, 0] / 2.0)
        }
        sum[0, 0] = 0.0
        return .circulant(vector: sum.elements) + .diagonal(vector: sum.fft1D().real.elements)
    }
}

extension ComplexMatrix where Scalar == Double {
    
    public static func dfrftEigenvalues(length: Int, order: Scalar) -> ComplexMatrix<Scalar> {
        let alpha = order * (.pi / 2)
        let indices = Matrix.dfrftIndices(length: length)
        return ComplexMatrix(real: (indices * -alpha).cos(), imaginary: (indices * -alpha).sin())
    }
    
}

extension Matrix where Scalar == Double {
    
    public static func dfrftIndices(length: Int) -> Matrix<Scalar> {
        return .init(shape: .row(length: length)) { row, column in
            return Scalar(column < length - 1 ? column : length - length % 2)
        }
    }
    
}

extension Matrix where Scalar == Double {
    
    public static func centralDifference(order: Int) -> Matrix {
        return .init(shape: .row(length: order + 1)) { row, column in
            let n = order
            let k = column
            return (k % 2 == 0 ? 1.0 : -1.0) * tgamma(Double(n + 1)) / (tgamma(Double(k + 1)) * tgamma(Double(n - k + 1)))
        }
    }
    
}

extension Matrix where Scalar == Double {
    
    internal static func oddEvenDecompositionMatrix(length: Int) -> Matrix<Scalar> {
        return .init(shape: .square(length: length)) { row, column in
            let diagonal = column - row
            let antidiagonal = column + row
            switch (diagonal, antidiagonal) {
            case (0, 0):
                return 1.0
            case (0, length):
                return 1.0
            case (0, _) where column + row < length:
                return 1.0 / Scalar.sqrt(2.0)
            case (0, _):
                return -1.0 / Scalar.sqrt(2.0)
            case (_, length):
                return 1.0 / Scalar.sqrt(2.0)
            default:
                return 0.0
            }
        }
    }
    
}
