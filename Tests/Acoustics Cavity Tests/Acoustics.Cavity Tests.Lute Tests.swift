import Acoustics_Cavity
import Acoustics_Plate
import Acoustics_Shell
import Mesh
import Testing

extension `Acoustics.Cavity Tests` {
    /// A lute's whole body from its geometry: a barred spruce top with its bridge, a bowl of ribs, and the air inside,
    /// the top and the bowl one shell joined at the rim, against Firth's (1977) measured admittance at the bridge.
    @Suite struct Lute {}
}

extension `Acoustics.Cavity Tests`.Lute {
    typealias Body = `Acoustics.Cavity Tests`.Body
    typealias Space3 = Point<3, Double>

    /// Gamut's Fantasia model's size, 43.2 cm long, 31.5 cm wide and 15.3 cm deep, as in `Body.lute`.
    static let length = 0.432
    static let centre = Point<2, Double>(x: -0.02, y: 0)
    static let depth = 0.153
    /// A point `centimetres` from the tail of Firth's 45 cm body, scaled to this one.
    static func along(_ centimetres: Double) -> Double { -0.18 + length * centimetres / 45 }
    /// The outline's half width at x, close to `Body.lute`'s.
    static func half(_ x: Double) -> Double {
        0.1575 * Swift.max(0, 1 - (x + 0.02) * (x + 0.02) / (0.24 * 0.24)).squareRoot()
    }

    static let spruce = Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter)
    static let walnut = Acoustics.Plate.Material(.blackWalnut, cut: .quarter)
    static let up = Vector<3, Double>(x: 0, y: 0, z: 1)

    /// How the bowl is taken.
    enum Bowl {
        /// Held still: only the top moves, clamped at the rim.
        case rigid
        /// A wall of its own beside the top, both clamped at the rim they share.
        case separate
        /// One shell with the top, joined at the rim, held only at the neck block.
        case joined
    }

    /// What a body is built from.
    struct Options {
        var top = 0.0025
        var ribs = 0.0025
        var bars = true
        var bridge = true
        var bowl = Bowl.joined
        var cells = 8
    }

    /// Firth's four bars at 10, 20, 30 and 40 cm, inside, with Schreiner's section from the Jauch lute: 4.2 mm wide,
    /// 15.8 mm high in the middle falling to 6 mm at the ribs, each cut into pieces of their mean height.
    static var bars: [Acoustics.Shell.Bar] {
        [10.0, 20, 30, 40].flatMap { centimetres -> [Acoustics.Shell.Bar] in
            let x = along(centimetres)
            let h = half(x)
            let count = 12
            return (0..<count).map { index in
                let (a, b) = (Double(index) / Double(count), Double(index + 1) / Double(count))
                let share = (abs(2 * (a + b) / 2 - 1))
                let height = 0.006 + 0.0098 * (1 - share * share)
                let reach = 1.15 * h
                return Acoustics.Shell.Bar(
                    from: Point(x: x, y: -reach + 2 * reach * a, z: 0),
                    to: Point(x: x, y: -reach + 2 * reach * b, z: 0),
                    width: Space.Physical.metres(0.0042),
                    height: Space.Physical.metres(height),
                    material: spruce,
                    side: Vector(x: 0, y: 0, z: -1)
                )
            }
        }
    }

    /// The bridge at Firth's 11.5 cm, glued on top: walnut, 130 mm long, 9 mm wide and 6 mm high, about 4.5 g.
    static var bridge: Acoustics.Shell.Bar {
        Acoustics.Shell.Bar(
            from: Point(x: along(11.5), y: -0.065, z: 0),
            to: Point(x: along(11.5), y: 0.065, z: 0),
            width: Space.Physical.metres(0.009),
            height: Space.Physical.metres(0.006),
            material: walnut,
            side: up
        )
    }

    /// A body: its coupled modes, the vertex at the treble end of the bridge, and the walls.
    struct Model {
        let body: Acoustics.Cavity.Coupling
        let bridge: Int
        let walls: [Acoustics.Cavity.Wall]
        let shells: [Acoustics.Shell.Modes]
    }

    static func model(_ options: Options) throws -> Model {
        let outline = Body.lute
        let volume = try Mesh<Space3>.Volume(
            under: outline,
            around: centre,
            depth: depth,
            cells: Cardinal(UInt(options.cells))
        )
        let cavity = Body.cavity(volume)
        let air = try Acoustics.Cavity.Modes(cavity, count: 12)
        let points = volume.vertices
        let rose = Point<2, Double>(x: along(35), y: 0)
        var top: [[Int]] = []
        var bowl: [[Int]] = []
        for face in volume.boundary.faces {
            let corners = face.vertices.map { Int(bitPattern: $0) }
            let p = corners.map { points[$0] }
            if p.allSatisfy({ abs($0.z) < 1e-12 }) {
                let cx = p.map { $0.x }.reduce(0, +) / 4
                let cy = p.map { $0.y }.reduce(0, +) / 4
                if (cx - rose.x) * (cx - rose.x) + (cy - rose.y) * (cy - rose.y) >= 0.04 * 0.04 { top.append(corners) }
            } else {
                bowl.append(corners)
            }
        }
        let rim = Set(top.flatMap { $0 }).intersection(Set(bowl.flatMap { $0 }))
        let neck = Set(points.indices.filter { points[$0].x > 0.232 })
        func shell(_ faces: [[Int]], ribs: Bool, top: Bool, held: Set<Int>) throws -> Acoustics.Shell {
            let mesh = try Mesh(
                vertices: points,
                faces: faces.map { face in try Polygon(vertices: face.map { Index(Ordinal(UInt($0))) }) }
            )
            var thicknesses: [Int: Measure<1, Space.Physical>.Value<Double>] = [:]
            var materials: [Int: Acoustics.Plate.Material] = [:]
            let topCount = top ? (ribs ? faces.count - bowl.count : faces.count) : 0
            for index in 0..<topCount {
                thicknesses[index] = Space.Physical.metres(options.top)
                materials[index] = spruce
            }
            var bars: [Acoustics.Shell.Bar] = []
            if top && options.bars { bars += Self.bars }
            if top && options.bridge { bars.append(Self.bridge) }
            let used = Set(faces.flatMap { $0 })
            let clamped = Acoustics.Shell.Hold.clamped
            return Acoustics.Shell(
                mesh: mesh,
                thickness: Space.Physical.metres(options.ribs),
                material: walnut,
                grain: Vector(x: 1, y: 0, z: 0),
                thicknesses: thicknesses,
                materials: materials,
                bars: bars,
                held: Dictionary(uniqueKeysWithValues: held.intersection(used).map { ($0, clamped) })
            )
        }
        func wall(_ modes: Acoustics.Shell.Modes, faces: [[Int]]) -> Acoustics.Cavity.Wall {
            Acoustics.Cavity.Wall(
                faces: faces,
                modes: modes.map { mode in
                    Acoustics.Cavity.Wall.Mode(
                        frequency: mode.frequency,
                        mass: mode.mass,
                        loss: 0.01,
                        displacement: mode.shape
                    )
                }
            )
        }
        var walls: [Acoustics.Cavity.Wall] = []
        var shells: [Acoustics.Shell.Modes] = []
        switch options.bowl {
        case .rigid:
            let alone = try shell(top, ribs: false, top: true, held: rim.union(neck))
            let modes = try Acoustics.Shell.Modes(alone, count: 20)
            shells = [modes]
            walls = [wall(modes, faces: top)]
        case .separate:
            let above = try shell(top, ribs: false, top: true, held: rim.union(neck))
            let below = try shell(bowl, ribs: true, top: false, held: rim.union(neck))
            let upper = try Acoustics.Shell.Modes(above, count: 20)
            let lower = try Acoustics.Shell.Modes(below, count: 30)
            shells = [upper, lower]
            walls = [wall(upper, faces: top), wall(lower, faces: bowl)]
        case .joined:
            let modes = try Acoustics.Shell.Modes(try shell(top + bowl, ribs: true, top: true, held: neck), count: 40)
            shells = [modes]
            walls = [wall(modes, faces: top + bowl)]
        }
        let opening = Acoustics.Cavity.Opening.circle(
            radius: Space.Physical.metres(0.04),
            thickness: Space.Physical.metres(options.top)
        )
        let body = try Acoustics.Cavity.Coupling(air, walls: walls, openings: [opening], count: 16)
        let target = (x: along(11.5), y: -0.05)
        let topVertices = Set(top.flatMap { $0 })
        let bridge = topVertices.min { a, b in
            let (p, q) = (points[a], points[b])
            return (p.x - target.x) * (p.x - target.x) + (p.y - target.y) * (p.y - target.y)
                < (q.x - target.x) * (q.x - target.x) + (q.y - target.y) * (q.y - target.y)
        }!
        return Model(body: body, bridge: bridge, walls: walls, shells: shells)
    }

    /// The peaks of the admittance's magnitude at the bridge, along z, from 60 to 800 Hz.
    static func peaks(_ model: Model) -> [Double] {
        let reach = model.body.reach(at: model.bridge, along: up)
        let frequencies = stride(from: 60.0, through: 800, by: 0.5).map { $0 }
        let magnitude = frequencies.map { Body.admittance(model.body, reach: reach, at: $0) }
        return (1..<(frequencies.count - 1)).filter {
            magnitude[$0] > magnitude[$0 - 1] && magnitude[$0] > magnitude[$0 + 1]
        }.map { frequencies[$0] }
    }

    /// How far each body mode moves the top along z and the bowl, at most, at its own scale.
    static func reach(_ model: Model) -> [(top: Double, bowl: Double)] {
        let shell = model.shells[0]
        let points = shell.shell.mesh.vertices
        return model.body.map { mode in
            var motion = [Vector<3, Double>](repeating: .zero, count: points.count)
            for (share, wallMode) in zip(mode.walls[0], shell) {
                for (index, v) in wallMode.shape.enumerated() { motion[index] += v.scaled(by: share) }
            }
            var top = 0.0
            var bowl = 0.0
            for (index, v) in motion.enumerated() {
                if abs(points[index].z) < 1e-12 {
                    top = Swift.max(top, abs(v.z))
                } else {
                    bowl = Swift.max(bowl, v.length)
                }
            }
            return (top, bowl)
        }
    }

    @Test
    func `The joined body's peaks at the bridge are near Firth's`() throws {
        // Firth (1977), the input admittance at the treble end of a lute's bridge (Rossing 2010, Fig. 3.26a, p. 41):
        // the air at 132 Hz, the top at 254, 304, 395 and 602 Hz. The lute is not Firth's.
        let model = try Self.model(Options(cells: 12))
        let peaks = Self.peaks(model)
        func nearest(_ f: Double) -> Double { peaks.min { abs($0 - f) < abs($1 - f) }! }
        #expect(abs(nearest(132) / 132 - 1) < 0.05, "\(peaks)")
        for measured in [304.0, 395, 602] {
            #expect(abs(nearest(measured) / measured - 1) < 0.1, "\(measured): \(peaks)")
        }
        // His 254 Hz has no mode of the top near it. The nearest peak is the body rocking on its neck block, the top
        // and the bowl moving alike and no air through the rose, about a fifth low.
        let rocking = nearest(254)
        #expect(abs(rocking / 254 - 1) < 0.25, "\(peaks)")
        let index = model.body.indices.min {
            abs(model.body[$0].frequency.underlying - rocking) < abs(model.body[$1].frequency.underlying - rocking)
        }!
        let moving = Self.reach(model)[index]
        #expect(moving.bowl > 0.5 * moving.top)
        #expect(abs(model.body[index].openings[0]) < 1e-6)
    }

    @Test
    func `Joining the bowl at the rim lowers the top's peaks, holding it does not`() throws {
        // With the rim clamped the top rings at its #345 frequencies, 400 Hz and up, whether the bowl is rigid or a
        // wall of its own; joined, the top's first peak above the air falls by a fifth.
        let rigid = Self.peaks(try Self.model(Options(bowl: .rigid)))
        let separate = Self.peaks(try Self.model(Options(bowl: .separate)))
        let joined = try Self.model(Options(bowl: .joined))
        let firstTop = { (peaks: [Double]) in peaks.first { $0 > 200 }! }
        #expect(firstTop(rigid) > 380, "\(rigid)")
        #expect(abs(firstTop(separate) / firstTop(rigid) - 1) < 0.02, "\(separate)")
        let top = Self.topPeak(joined)
        #expect(top < 0.85 * firstTop(rigid), "\(top)")
    }

    /// The first peak above 200 Hz whose mode moves the top more than twice as far as the bowl.
    static func topPeak(_ model: Model) -> Double {
        let reach = Self.reach(model)
        let peaks = Self.peaks(model)
        return peaks.first { peak in
            guard peak > 200 else { return false }
            let index = model.body.indices.min {
                abs(model.body[$0].frequency.underlying - peak) < abs(model.body[$1].frequency.underlying - peak)
            }!
            return reach[index].top > 2 * reach[index].bowl
        }!
    }
}
