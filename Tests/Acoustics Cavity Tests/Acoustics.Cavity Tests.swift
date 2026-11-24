import Acoustics_Cavity
import Testing

@Suite struct `Acoustics.Cavity Tests` {}

extension `Acoustics.Cavity Tests` {
    /// Structured meshes of the solids the closed forms are known for.
    enum Solid {
        /// The box [0, x] × [0, y] × [0, z] cut into `cells` hexahedra along each axis, or each hexahedron cut
        /// into six tetrahedra around its diagonal from corner 0 to corner 6.
        static func box(
            _ x: Double,
            _ y: Double,
            _ z: Double,
            cells: (Int, Int, Int),
            tetrahedra: Bool = false
        ) -> Acoustics.Cavity {
            let (nx, ny, nz) = cells
            var vertices: [Point<3, Double>] = []
            for k in 0...nz {
                for j in 0...ny {
                    for i in 0...nx {
                        let u = Double(i) / Double(nx)
                        let v = Double(j) / Double(ny)
                        let w = Double(k) / Double(nz)
                        vertices.append(Point(x: x * u, y: y * v, z: z * w))
                    }
                }
            }
            func index(_ i: Int, _ j: Int, _ k: Int) -> Int { (k * (ny + 1) + j) * (nx + 1) + i }
            var result: [[Int]] = []
            for k in 0..<nz {
                for j in 0..<ny {
                    for i in 0..<nx {
                        let c = [
                            index(i, j, k), index(i + 1, j, k), index(i + 1, j + 1, k), index(i, j + 1, k),
                            index(i, j, k + 1), index(i + 1, j, k + 1), index(i + 1, j + 1, k + 1),
                            index(i, j + 1, k + 1),
                        ]
                        if tetrahedra {
                            for path in [[1, 2], [2, 3], [3, 7], [7, 4], [4, 5], [5, 1]] {
                                result.append([c[0], c[path[0]], c[path[1]], c[6]])
                            }
                        } else {
                            result.append(c)
                        }
                    }
                }
            }
            return Acoustics.Cavity(vertices: vertices, cells: result)
        }

        /// A disk's section: a core square of `core` × `core` cells (each side 2 `core` points apart, so `core` is
        /// half its cells across), then `rings` rings carrying its outline to a circle. Ring r of `rings` lies at
        /// `radius(r)`, ring 0 being the core's outline; a radius of nil blends the square into the circle of the
        /// first ring that has one.
        struct Section {
            let core: Int
            let rings: Int

            /// The half-width of the grid in points.
            var half: Int { core + rings }

            /// The ring of grid point (i, j), 0 inside the core.
            func ring(_ i: Int, _ j: Int) -> Int { Swift.max(Swift.max(Swift.abs(i), Swift.abs(j)) - core, 0) }

            /// The point (i, j) of the section, with the core's half-side `side` and ring radii from `radius`.
            func point(_ i: Int, _ j: Int, side: Double, blend: Int, radius: (Int) -> Double) -> (Double, Double) {
                let level = Swift.max(Swift.abs(i), Swift.abs(j))
                let u = Double(i) / Double(Swift.max(level, 1))
                let v = Double(j) / Double(Swift.max(level, 1))
                let length = (u * u + v * v).squareRoot()
                if level <= core { return (side * Double(i) / Double(core), side * Double(j) / Double(core)) }
                let r = level - core
                if r <= blend {
                    let t = Double(r) / Double(blend)
                    let circle = radius(blend)
                    return ((1 - t) * side * u + t * circle * u / length, (1 - t) * side * v + t * circle * v / length)
                }
                return (radius(r) * u / length, radius(r) * v / length)
            }

            /// The section's quads, counterclockwise, as (i, j) corners, keeping those whose rings are all at most
            /// `within`.
            func quads(within: Int) -> [[(Int, Int)]] {
                var result: [[(Int, Int)]] = []
                for j in -half..<half {
                    for i in -half..<half {
                        let quad = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                        if quad.allSatisfy({ ring($0.0, $0.1) <= within }) { result.append(quad) }
                    }
                }
                return result
            }
        }

        /// Stacks of a section along z: each slab between two levels keeps the quads whose rings are at most its
        /// `within`, and each vertex takes its section point from the radii of the first slab that uses it.
        static func stack(
            _ section: Section,
            levels: [Double],
            within: [Int],
            point: (Int, Int, Int) -> (Double, Double),
            open: (Int, Int, Int) -> Bool = { _, _, _ in false }
        ) -> Acoustics.Cavity {
            var numbers: [Int: Int] = [:]
            var vertices: [Point<3, Double>] = []
            var opened: Set<Int> = []
            let width = 2 * section.half + 1
            func number(_ i: Int, _ j: Int, _ k: Int) -> Int {
                let key = (k * width + (j + section.half)) * width + (i + section.half)
                if let found = numbers[key] { return found }
                let (x, y) = point(i, j, k)
                vertices.append(Point(x: x, y: y, z: levels[k]))
                numbers[key] = vertices.count - 1
                if open(i, j, k) { opened.insert(vertices.count - 1) }
                return vertices.count - 1
            }
            var cells: [[Int]] = []
            for k in 0..<(levels.count - 1) {
                for quad in section.quads(within: within[k]) {
                    cells.append(quad.map { number($0.0, $0.1, k) } + quad.map { number($0.0, $0.1, k + 1) })
                }
            }
            return Acoustics.Cavity(vertices: vertices, cells: cells, open: opened)
        }

        /// The cylinder of `radius` about the z axis from z = 0 to `height`.
        static func cylinder(radius: Double, height: Double, core: Int, rings: Int, layers: Int) -> Acoustics.Cavity {
            let section = Section(core: core, rings: rings)
            return stack(
                section,
                levels: (0...layers).map { height * Double($0) / Double(layers) },
                within: Array(repeating: rings, count: layers),
                point: { i, j, _ in section.point(i, j, side: 0.5 * radius, blend: rings) { _ in radius } }
            )
        }

        /// The ball of `radius` about the origin: a core cube of half-side 0.35 `radius` with `core` cells along each
        /// half-axis, and `shells` layers carrying each of its faces along rays onto the sphere.
        static func ball(radius: Double, core: Int, shells: Int) -> Acoustics.Cavity {
            let half = core + shells
            let width = 2 * half + 1
            var vertices: [Point<3, Double>] = []
            for k in -half...half {
                for j in -half...half {
                    for i in -half...half {
                        let level = Swift.max(Swift.abs(i), Swift.abs(j), Swift.abs(k))
                        let side = 0.35 * radius
                        if level <= core {
                            let s = side / Double(core)
                            vertices.append(Point(x: s * Double(i), y: s * Double(j), z: s * Double(k)))
                            continue
                        }
                        let u = [Double(i), Double(j), Double(k)].map { $0 / Double(level) }
                        let length = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).squareRoot()
                        let t = Double(level - core) / Double(shells)
                        let p = u.map { (1 - t) * side * $0 + t * radius * $0 / length }
                        vertices.append(Point(x: p[0], y: p[1], z: p[2]))
                    }
                }
            }
            func index(_ i: Int, _ j: Int, _ k: Int) -> Int {
                ((k + half) * width + (j + half)) * width + (i + half)
            }
            var cells: [[Int]] = []
            for k in -half..<half {
                for j in -half..<half {
                    for i in -half..<half {
                        cells.append([
                            index(i, j, k), index(i + 1, j, k), index(i + 1, j + 1, k), index(i, j + 1, k),
                            index(i, j, k + 1), index(i + 1, j, k + 1), index(i + 1, j + 1, k + 1),
                            index(i, j + 1, k + 1),
                        ])
                    }
                }
            }
            return Acoustics.Cavity(vertices: vertices, cells: cells)
        }
    }

    /// A Helmholtz resonator about the z axis: a cylindrical body of `body` radius and depth below z = 0, a
    /// neck of `neck` radius and length above it, and the air above the neck's flange out to a cylinder of
    /// `outside` radius and height, released there.
    static func resonator(
        body: (radius: Double, depth: Double),
        neck: (radius: Double, length: Double),
        outside: (radius: Double, height: Double),
        fineness: Int = 2
    ) -> Acoustics.Cavity {
        let blend = (3 * fineness + 1) / 2
        let rings = blend + 3 * fineness
        let section = Solid.Section(core: (3 * fineness + 1) / 2, rings: rings)
        /// The growth g of steps that start at `first` and add up to `span` in `count` steps, by halving.
        func growth(first: Double, span: Double, count: Int) -> Double {
            var low = 1.0
            var high = 10.0
            for _ in 0..<100 {
                let middle = (low + high) / 2
                var sum = 0.0
                var step = first
                for _ in 0..<count {
                    sum += step
                    step *= middle
                }
                if sum < span { low = middle } else { high = middle }
            }
            return low
        }
        func graded(_ length: Double, _ count: Int, _ exponent: Int) -> [Double] {
            (0...count).map { k in
                var share = 1.0
                for _ in 0..<exponent { share *= Double(k) / Double(count) }
                return length * share
            }
        }
        let a = neck.radius
        let below = graded(body.depth, 4 * fineness, 2).reversed().map { -$0 }
        let along = graded(neck.length, 4 * fineness, 1).dropFirst()
        let above = graded(outside.height, 5 * fineness, 2).dropFirst().map { neck.length + $0 }
        let levels = below + along + above
        let within = levels.indices.dropLast().map { k in
            levels[k] >= 0 && levels[k + 1] <= neck.length ? blend : rings
        }
        return Solid.stack(
            section,
            levels: levels,
            within: within,
            point: { i, j, k in
                // Beyond the neck the rings start at the blend's own spacing and grow geometrically.
                let outer = levels[k] <= 0 ? body.radius : outside.radius
                let first = 0.5 * a / Double(blend)
                let g = growth(first: first, span: outer - a, count: rings - blend)
                return section.point(i, j, side: 0.5 * a, blend: blend) { r in
                    var radius = a
                    var step = first
                    if r > blend {
                        for _ in 0..<(r - blend) {
                            radius += step
                            step *= g
                        }
                    }
                    return radius
                }
            },
            open: { i, j, k in
                k == levels.count - 1 || (levels[k] >= neck.length && section.ring(i, j) == rings)
            }
        )
    }

    static let air = Acoustics.Air.room
    static var c: Double { air.speed.underlying }

    /// The relative difference.
    static func error(_ found: Double, _ exact: Double) -> Double { Swift.abs(found - exact) / exact }
}
