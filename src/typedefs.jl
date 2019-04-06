const SymAnyDict = Dict{Symbol, Any}

const Point{N, F} = SVector{N, F}
const Point3{F} = SVector{3, F}
const Vec3{F}   = SVector{3, F}

const Mat3{F}   = SMatrix{3, 3, F, 9}
const Mat4{F}   = SMatrix{4, 4, F, 16}

export Point, Point3, Vec3, Mat3, Mat4

