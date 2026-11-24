public import Eigen

extension Acoustics.Cavity {
    /// Why a cavity's modes could not be found.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// The cell at this index has neither eight nor four corners, names a vertex that does not exist, or turns
        /// inside out somewhere.
        case cell(Int)
        /// The eigen solver failed.
        case eigen(Eigen.Error)
        /// More modes were asked for than the cavity has.
        case count
        /// The wall at this index covers a face that is neither a triangle nor a quadrilateral of the cavity's
        /// vertices, or has a mode whose displacement is not given at every vertex.
        case wall(Int)
        /// An opening was given, but none of the air's modes is the cavity's uniform pressure: the cavity is not
        /// closed, or its zero mode was not computed.
        case opening
    }
}
