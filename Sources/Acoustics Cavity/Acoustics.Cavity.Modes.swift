public import Acoustics
public import Eigen
internal import Galerkin

extension Acoustics.Cavity {
    /// A body of air's lowest acoustic modes, computed once from its shape and kept as a value.
    public struct Modes {
        /// The cavity they belong to.
        public let cavity: Acoustics.Cavity
        /// The air that fills it.
        public let air: Acoustics.Air
        /// The modes, lowest first; a repeated frequency appears once for each mode that has it.
        public let elements: [Mode]

        /// The `count` lowest modes of `cavity` filled with `air`.
        ///
        /// A cavity with no open vertex has a first mode at zero frequency to rounding: a pressure the same
        /// everywhere.
        public init(_ cavity: Acoustics.Cavity, count: Int, air: Acoustics.Air = .room) throws(Acoustics.Cavity.Error) {
            let system = try Acoustics.Cavity.System(cavity)
            guard count >= 0, count <= system.order else { throw .count }
            let pairs: [Eigen.Pair<Double>]
            do throws(Eigen.Error) {
                let pencil = try Eigen.Pencil(system.stiffness, system.mass)
                pairs = try pencil.lowest(Cardinal(UInt(count)))
            } catch {
                throw .eigen(error)
            }
            self.cavity = cavity
            self.air = air
            self.elements = pairs.map { system.mode($0, speed: air.speed.underlying) }
        }
    }
}

extension Acoustics.Cavity.Modes: Hashable {}
extension Acoustics.Cavity.Modes: Sendable {}

extension Acoustics.Cavity.Modes: RandomAccessCollection {
    public var startIndex: Int { elements.startIndex }
    public var endIndex: Int { elements.endIndex }
    public subscript(position: Int) -> Acoustics.Cavity.Mode { elements[position] }
}

extension Acoustics.Cavity {
    /// The weak form of the Helmholtz equation over the cavity's unknown pressures: ∫ ∇N·∇N dV p = k² ∫ N N dV p,
    /// with k = ω / c. A rigid wall is its natural boundary condition; an open vertex is held at zero.
    internal struct System {
        let stiffness: Sparse.Symmetric<Double>
        let mass: Sparse.Symmetric<Double>
        /// For each vertex, the unknown its pressure is, or nil where it is open.
        let pressure: [Int?]

        var order: Int { stiffness.order }
    }
}

extension Acoustics.Cavity.System {
    /// A reference cell's shape functions at its Gauss points, which every cell of that kind shares.
    struct Rule {
        let shapes: [Galerkin.Shape<Double>]
        let weights: [Double]
        /// The corners in the other orientation.
        let turned: [Int]
    }

    init(_ cavity: Acoustics.Cavity) throws(Acoustics.Cavity.Error) {
        let points = cavity.vertices
        // Degree 3 on the cube and 2 on the tetrahedron integrate the mass matrix exactly on an affine cell.
        let cube = Galerkin.Reference.Hexahedron.quadrature(degree: 3, Double.self)
        let simplex = Galerkin.Reference.Tetrahedron.quadrature(degree: 2, Double.self)
        let hexahedron = Rule(
            shapes: cube.points.map { Galerkin.Reference.Hexahedron.shape(at: $0) },
            weights: cube.weights,
            turned: [4, 5, 6, 7, 0, 1, 2, 3]
        )
        let tetrahedron = Rule(
            shapes: simplex.points.map { Galerkin.Reference.Tetrahedron.shape(at: $0) },
            weights: simplex.weights,
            turned: [0, 2, 1, 3]
        )

        let unknowns = Galerkin.Unknowns(
            vertices: points.count,
            components: 1,
            held: Dictionary(uniqueKeysWithValues: cavity.open.map { ($0, [0]) })
        )
        var stiffness = Galerkin.Assembly<Double>(unknowns)
        var mass = stiffness
        for (index, cell) in cavity.cells.enumerated() {
            let rule: Rule
            switch cell.count {
            case 8: rule = hexahedron
            case 4: rule = tetrahedron
            default: throw .cell(index)
            }
            guard cell.allSatisfy({ points.indices.contains($0) }) else { throw .cell(index) }
            var corners = cell
            var jacobians = rule.shapes.map { Galerkin.Jacobian.Solid($0, corners: corners.map { points[$0] }) }
            if jacobians.allSatisfy({ $0.determinant < 0 }) {
                corners = rule.turned.map { cell[$0] }
                jacobians = rule.shapes.map { Galerkin.Jacobian.Solid($0, corners: corners.map { points[$0] }) }
            }
            guard jacobians.allSatisfy({ $0.determinant > 0 }) else { throw .cell(index) }
            let n = corners.count
            var k = [Double](repeating: 0, count: n * n)
            var m = k
            for (point, (shape, jacobian)) in zip(rule.weights, zip(rule.shapes, jacobians)) {
                let volume = point * jacobian.determinant
                let gradient = jacobian.gradients(shape)
                for a in 0..<n {
                    for b in a..<n {
                        let dot =
                            gradient.x[a] * gradient.x[b] + gradient.y[a] * gradient.y[b]
                            + gradient.z[a] * gradient.z[b]
                        k[n * a + b] += volume * dot
                        m[n * a + b] += volume * shape.values[a] * shape.values[b]
                    }
                }
            }
            for a in 0..<n {
                for b in 0..<a {
                    k[n * a + b] = k[n * b + a]
                    m[n * a + b] = m[n * b + a]
                }
            }
            stiffness.add(k, at: corners)
            mass.add(m, at: corners)
        }
        self.stiffness = stiffness.matrix
        self.mass = mass.matrix
        self.pressure = points.indices.map { unknowns[$0, 0] }
    }

    /// The mode of an eigenpair: its frequency, its shape scaled to a largest pressure of one, and its modal volume
    /// at that scale.
    func mode(_ pair: Eigen.Pair<Double>, speed: Double) -> Acoustics.Cavity.Mode {
        var shape = pressure.map { $0.map { pair.vector[$0] } ?? 0 }
        var peak = 0.0
        for value in shape where Swift.abs(value) > Swift.abs(peak) { peak = value }
        if peak != 0 {
            for index in shape.indices { shape[index] /= peak }
        }
        // The eigenvalue is k², and the vector is normalized so that ∫ p² dV = 1, so at a largest pressure of one
        // the modal volume is 1 / peak².
        let frequency = speed * Swift.max(pair.value, 0).squareRoot() / (2 * Double.pi)
        return Acoustics.Cavity.Mode(
            frequency: Frequency.hertz(frequency),
            shape: shape,
            volume: Volume<Space.Physical>.amount(peak != 0 ? 1 / (peak * peak) : 0)
        )
    }
}
