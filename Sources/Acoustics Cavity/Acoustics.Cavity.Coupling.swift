public import Acoustics
public import Eigen
internal import Galerkin

extension Acoustics.Cavity {
    /// A cavity's air coupled to the walls that move and the openings it breathes through: the modes of a whole
    /// body, such as a guitar's or a lute's, computed once from the modes of its parts and kept as a value.
    ///
    /// The air's modes and each wall's modes in vacuo are the basis. Pressure on a wall drives its modes, and a
    /// wall's motion drives the air's: mode j's amplitude a, the pressure a × shape, obeys
    /// ä + ωⱼ² a = −(ρc² / Λⱼ) ∫ shape ü dS. The cavity's zero mode, a pressure the same everywhere, is a spring
    /// against the volume the walls and the openings push out, (ρc² / V)(Σ swept volume)²; each other mode adds the
    /// air's mass to the walls and its own resonance. The coupled problem is symmetric: K x = ω² M x with
    ///
    ///     K = [Kʷ + K⁰   0 ]     M = [Mʷ + D Dᵀ   D]
    ///         [0         Ω²]         [Dᵀ          I]
    ///
    /// over the walls' and openings' amplitudes and the other air modes' scaled amplitudes, where
    /// Dₖⱼ = ∫ φₖ ψⱼ dS √(ρc² / Λⱼ) / ωⱼ.
    public struct Coupling {
        /// The air's modes the coupling is built on.
        public let air: Acoustics.Cavity.Modes
        /// The walls that move.
        public let walls: [Wall]
        /// The openings.
        public let openings: [Opening]
        /// The coupled modes, lowest first.
        public let elements: [Mode]

        /// The `count` lowest modes of `air` coupled to `walls` and `openings`.
        ///
        /// The more air and wall modes go in, the closer the low coupled modes come to the body's; the highest
        /// coupled modes are only as good as the basis allows.
        public init(
            _ air: Acoustics.Cavity.Modes,
            walls: [Wall],
            openings: [Opening] = [],
            count: Int
        ) throws(Acoustics.Cavity.Error) {
            let system = try System(air, walls: walls, openings: openings)
            guard count >= 0, count <= system.order else { throw .count }
            let pairs: [Eigen.Pair<Double>]
            do throws(Eigen.Error) {
                let pencil = try Eigen.Pencil(system.stiffness, system.mass)
                pairs = try pencil.lowest(Cardinal(UInt(count)))
            } catch {
                throw .eigen(error)
            }
            self.air = air
            self.walls = walls
            self.openings = openings
            self.elements = pairs.map { system.mode($0) }
        }
    }
}

extension Acoustics.Cavity.Coupling: Hashable {}
extension Acoustics.Cavity.Coupling: Sendable {}

extension Acoustics.Cavity.Coupling: RandomAccessCollection {
    public var startIndex: Int { elements.startIndex }
    public var endIndex: Int { elements.endIndex }
    public subscript(position: Int) -> Mode { elements[position] }
}

extension Acoustics.Cavity.Coupling {
    /// How far each mode moves `vertex` along `direction`, a unit vector, per unit of its amplitude: the sum of its
    /// walls' modes' motions there, weighted by their shares. A vertex no wall covers does not move.
    ///
    /// Struck there along `direction`, each mode is a resonator of mass `mass / r²`.
    public func reach(at vertex: Int, along direction: Vector<3, Double>) -> [Double] {
        elements.map { mode in
            var sum = 0.0
            for (wall, shares) in zip(walls, mode.walls) {
                for (share, wallMode) in zip(shares, wall.modes) where vertex < wallMode.displacement.count {
                    sum += share * wallMode.displacement[vertex].dot(direction)
                }
            }
            return sum
        }
    }
}

extension Acoustics.Cavity.Coupling {
    /// The coupled problem's stiffness and mass, and what turns an eigenvector back into a mode.
    internal struct System {
        let stiffness: Sparse.Symmetric<Double>
        let mass: Sparse.Symmetric<Double>
        /// Per wall mode, in order over the walls: its stiffness ω² m and loss factor.
        let wallStiffness: [Double]
        let wallLoss: [Double]
        /// How many amplitudes each wall has.
        let wallCounts: [Int]
        let openingCount: Int
        /// The volume each wall mode and each opening pushes out per unit of its amplitude.
        let swept: [Double]
        let density: Double
        let speed: Double

        var order: Int { stiffness.order }
    }
}

extension Acoustics.Cavity.Coupling.System {
    init(
        _ air: Acoustics.Cavity.Modes,
        walls: [Acoustics.Cavity.Wall],
        openings: [Acoustics.Cavity.Opening]
    ) throws(Acoustics.Cavity.Error) {
        let cavity = air.cavity
        let points = cavity.vertices
        let density = air.air.density.underlying
        let speed = air.air.speed.underlying
        let modes = air.elements
        // A mode is the cavity's uniform pressure when its wavelength is beyond a thousand times the cavity's
        // size, which a closed cavity's zero mode is to rounding.
        var size = 0.0
        if let first = points.first {
            var (low, high) = (first.coordinates, first.coordinates)
            for point in points {
                low = Vector(x: Swift.min(low.x, point.x), y: Swift.min(low.y, point.y), z: Swift.min(low.z, point.z))
                high = Vector(
                    x: Swift.max(high.x, point.x),
                    y: Swift.max(high.y, point.y),
                    z: Swift.max(high.z, point.z)
                )
            }
            size = (high - low).length
        }
        let isStatic = modes.map { 2 * Double.pi * $0.frequency.underlying * size / speed < 1e-3 }
        guard openings.isEmpty || isStatic.contains(true) else { throw .opening }

        // ∫ N N n dS over every wall's faces, n the unit normal out of the air, applied to each air mode's shape and
        // to one: vectors per vertex of the cavity.
        let square = Galerkin.Reference.Quadrilateral.quadrature(degree: 3, Double.self)
        let triangle = Galerkin.Reference.Triangle.quadrature(degree: 2, Double.self)
        let squareShapes = square.points.map { Galerkin.Reference.Quadrilateral.shape(at: $0) }
        let triangleShapes = triangle.points.map { Galerkin.Reference.Triangle.shape(at: $0) }

        var wallStiffness: [Double] = []
        var wallLoss: [Double] = []
        var wallMass: [Double] = []
        var swept: [Double] = []
        // coupling[k][j] = ∫ φₖ ψⱼ dS
        var coupling: [[Double]] = []
        for (index, wall) in walls.enumerated() {
            var weighted = [[Vector<3, Double>]](
                repeating: [Vector<3, Double>](repeating: .zero, count: points.count),
                count: modes.count
            )
            var area = [Vector<3, Double>](repeating: .zero, count: points.count)
            for face in wall.faces {
                guard face.count == 4 || face.count == 3, face.allSatisfy({ points.indices.contains($0) }) else {
                    throw .wall(index)
                }
                let shapes = face.count == 4 ? squareShapes : triangleShapes
                let weights = face.count == 4 ? square.weights : triangle.weights
                let corners = face.map { points[$0] }
                for (shape, weight) in zip(shapes, weights) {
                    let dS = Galerkin.Jacobian.Surface(shape, corners: corners).normal.scaled(by: weight)
                    for (a, vertex) in face.enumerated() {
                        let na = dS.scaled(by: shape.values[a])
                        area[vertex] += na
                        for (j, mode) in modes.enumerated() {
                            var pressure = 0.0
                            for (b, other) in face.enumerated() { pressure += shape.values[b] * mode.shape[other] }
                            weighted[j][vertex] += na.scaled(by: pressure)
                        }
                    }
                }
            }
            for mode in wall.modes {
                guard mode.displacement.count == points.count else { throw .wall(index) }
                let omega = 2 * Double.pi * mode.frequency.underlying
                wallMass.append(mode.mass.underlying)
                wallStiffness.append(omega * omega * mode.mass.underlying)
                wallLoss.append(mode.loss)
                var volume = 0.0
                for (d, s) in zip(mode.displacement, area) { volume += d.dot(s) }
                swept.append(volume)
                coupling.append(
                    weighted.map { row in
                        var sum = 0.0
                        for (d, w) in zip(mode.displacement, row) { sum += d.dot(w) }
                        return sum
                    }
                )
            }
        }
        let walled = wallMass.count
        for opening in openings {
            swept.append(opening.area.underlying)
            wallMass.append(density * opening.area.underlying * opening.length.underlying)
        }
        let structural = wallMass.count
        let dynamic = modes.indices.filter { !isStatic[$0] }
        let order = structural + dynamic.count

        var stiffness: [Sparse.Symmetric<Double>.Entry] = []
        var mass: [Sparse.Symmetric<Double>.Entry] = []
        func entry(_ row: Int, _ column: Int, _ value: Double) -> Sparse.Symmetric<Double>.Entry {
            Sparse.Symmetric<Double>.Entry(Ordinal(UInt(row)), Ordinal(UInt(column)), value)
        }
        // Kʷ + K⁰: each wall mode's own stiffness, and the uniform pressure against the volume pushed out.
        var uniform = [[Double]](repeating: [Double](repeating: 0, count: structural), count: structural)
        for j in modes.indices where isStatic[j] {
            let b = density * speed * speed / modes[j].volume.underlying
            // ∫ φ ψ₀ dS for the walls; S ψ₀ for an opening, ψ₀ being the same everywhere.
            let level = modes[j].shape.first ?? 1
            let g = (0..<structural).map { k in k < walled ? coupling[k][j] : swept[k] * level }
            for k in 0..<structural {
                for l in 0..<structural { uniform[k][l] += b * g[k] * g[l] }
            }
        }
        // D, and D Dᵀ on the walls' block.
        var d = [[Double]](repeating: [Double](repeating: 0, count: dynamic.count), count: walled)
        for (column, j) in dynamic.enumerated() {
            let omega = 2 * Double.pi * modes[j].frequency.underlying
            let scale = (density * speed * speed / modes[j].volume.underlying).squareRoot() / omega
            for k in 0..<walled { d[k][column] = coupling[k][j] * scale }
        }
        for k in 0..<structural {
            for l in k..<structural {
                let own = k == l && k < walled ? wallStiffness[k] : 0
                let value = own + uniform[k][l]
                if value != 0 { stiffness.append(entry(k, l, value)) }
                var m = k == l ? wallMass[k] : 0
                if k < walled && l < walled {
                    for column in dynamic.indices { m += d[k][column] * d[l][column] }
                }
                if m != 0 { mass.append(entry(k, l, m)) }
            }
            if k < walled {
                for column in dynamic.indices where d[k][column] != 0 {
                    mass.append(entry(k, structural + column, d[k][column]))
                }
            }
        }
        for (column, j) in dynamic.enumerated() {
            let omega = 2 * Double.pi * modes[j].frequency.underlying
            stiffness.append(entry(structural + column, structural + column, omega * omega))
            mass.append(entry(structural + column, structural + column, 1))
        }
        do throws(Sparse.Error) {
            self.stiffness = try Sparse.Symmetric(order: Cardinal(UInt(order)), entries: stiffness)
            self.mass = try Sparse.Symmetric(order: Cardinal(UInt(order)), entries: mass)
        } catch {
            throw .count
        }
        self.wallStiffness = wallStiffness
        self.wallLoss = wallLoss
        self.wallCounts = walls.map { $0.modes.count }
        self.openingCount = openings.count
        self.swept = swept
        self.density = density
        self.speed = speed
    }

    /// The mode of an eigenpair: its frequency, its shares scaled so that the largest amplitude of a wall mode or
    /// an opening is one, its mass at that scale, the volume it pushes out and its loss.
    func mode(_ pair: Eigen.Pair<Double>) -> Acoustics.Cavity.Coupling.Mode {
        let structural = swept.count
        var amplitudes = Array(pair.vector.prefix(structural))
        var peak = 0.0
        for value in amplitudes where Swift.abs(value) > Swift.abs(peak) { peak = value }
        if peak != 0 {
            for index in amplitudes.indices { amplitudes[index] /= peak }
        }
        let omega = Swift.max(pair.value, 0).squareRoot()
        // The vector is M-normalized, so at a largest amplitude of one the modal mass is 1 / peak².
        let mass = peak != 0 ? 1 / (peak * peak) : 0
        var volume = 0.0
        for (share, swept) in zip(amplitudes, swept) { volume += share * swept }
        // The walls' own losses, each weighted by the share of the strain energy it holds: Σ η k q² / (ω² m).
        var own = 0.0
        for (index, (k, eta)) in zip(wallStiffness, wallLoss).enumerated() {
            own += eta * k * amplitudes[index] * amplitudes[index]
        }
        // A small body radiates as a simple source of volume velocity ω V into free space, whose power over the
        // energy per radian of cycle is ρ ω V² / (4π m c).
        let loss: Double
        if mass > 0 && omega > 0 {
            loss = own / (omega * omega * mass) + density * omega * volume * volume / (4 * Double.pi * mass * speed)
        } else {
            loss = 0
        }
        var walls: [[Double]] = []
        var start = 0
        for count in wallCounts {
            walls.append(Array(amplitudes[start..<(start + count)]))
            start += count
        }
        return Acoustics.Cavity.Coupling.Mode(
            frequency: Frequency.hertz(omega / (2 * Double.pi)),
            mass: Mass.kilograms(mass),
            loss: loss,
            volume: Volume<Space.Physical>.amount(volume),
            walls: walls,
            openings: Array(amplitudes[start..<(start + openingCount)])
        )
    }
}
