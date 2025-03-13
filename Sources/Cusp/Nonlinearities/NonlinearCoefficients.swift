//
//  NonlinearCoefficients.swift
//  Cusp
//
//  Created by Daniel Clelland on 13/03/2025.
//

import Foundation
import Numerics
import Plinth

public struct NonlinearCoefficients<Scalar> where Scalar: Real {
    
    public var chi1: Scalar
    public var chi2: Scalar
    public var chi3: Scalar
    public var chi4: Scalar
    public var chi5: Scalar
    
    public init(chi1: Scalar = .zero, chi2: Scalar = .zero, chi3: Scalar = .zero, chi4: Scalar = .zero, chi5: Scalar = .zero) {
        self.chi1 = chi1
        self.chi2 = chi2
        self.chi3 = chi3
        self.chi4 = chi4
        self.chi5 = chi5
    }
    
    public subscript(_ order: Int) -> Scalar {
        get {
            switch order {
                case 1: return chi1
                case 2: return chi2
                case 3: return chi3
                case 4: return chi4
                case 5: return chi5
                default: fatalError()
            }
        }
        set {
            switch order {
                case 1: chi1 = newValue
                case 2: chi2 = newValue
                case 3: chi3 = newValue
                case 4: chi4 = newValue
                case 5: chi5 = newValue
                default: fatalError()
            }
        }
    }
    
}

extension ComplexMatrix where Scalar == Float {
    
    public func addNonlinearPhase(using coefficients: NonlinearCoefficients<Scalar>) -> ComplexMatrix {
        var phaseShift = Matrix.zeros(shape: shape)
        
        let amplitude = absolute()
        let intensity = squareMagnitudes()

        for order in 1...5 {
            if coefficients[order] != 0 {
                switch order {
                    case 1: phaseShift += coefficients[order]
                    case 2: phaseShift += coefficients[order] * amplitude
                    case 3: phaseShift += coefficients[order] * intensity
                    case 4: phaseShift += coefficients[order] * (amplitude * intensity)
                    case 5: phaseShift += coefficients[order] * (intensity * intensity)
                    default: break
                }
            }
        }

        return self * ComplexMatrix(imaginary: phaseShift).exp()
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func addNonlinearPhase(using coefficients: NonlinearCoefficients<Scalar>) -> ComplexMatrix {
        var phaseShift = Matrix.zeros(shape: shape)
        
        let amplitude = absolute()
        let intensity = squareMagnitudes()

        for order in 1...5 {
            if coefficients[order] != 0 {
                switch order {
                    case 1: phaseShift += coefficients[order]
                    case 2: phaseShift += coefficients[order] * amplitude
                    case 3: phaseShift += coefficients[order] * intensity
                    case 4: phaseShift += coefficients[order] * (amplitude * intensity)
                    case 5: phaseShift += coefficients[order] * (intensity * intensity)
                    default: break
                }
            }
        }

        return self * ComplexMatrix(imaginary: phaseShift).exp()
    }
    
}
