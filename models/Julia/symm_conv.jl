using NNlib: conv, ∇conv_data, depthwiseconv, output_size
using Flux: conv_reshape_bias, conv_dims, @functor

# The data structure
# The parameters have to be views into the weight array
# since many array elements share their values
struct Conv_Stencil{N,M,F,A,V}
  σ::F
  weight::A
  bias::V
  stride::NTuple{N,Int}
  pad::NTuple{M,Int}
  dilation::NTuple{N,Int}
  groups::Int
end

# The constructor
function Conv(k::NTuple{N,Integer}, ch::Pair{<:Integer,<:Integer}, σ = identity;
  init = glorot_uniform, stride = 1, pad = 0, dilation = 1, groups = 1,
  bias = true) where N

weight = convfilter(k, ch; init, groups)
Conv(weight, bias, σ; stride, pad, dilation, groups)
end

@functor Conv_Stencil

function (c::Conv)(x::AbstractArray)
  σ = NNlib.fast_act(c.σ, x)
  weight = gen_weight_array(c.weight)
  cdims = conv_dims(c, x)
  σ.(conv(x, weight, cdims) .+ conv_reshape_bias(c))
end

# Helper function
function gen_weight_array(x)

end

# Pretty printing from here on
function Base.show(io::IO, l::Conv_Stencil)
  print(io, "Conv(", size(l.weight)[1:ndims(l.weight)-2])
  print(io, ", ", _channels_in(l), " => ", _channels_out(l))
  _print_conv_opt(io, l)
  print(io, ")")
end

function _print_conv_opt(io::IO, l)
  l.σ == identity || print(io, ", ", l.σ)
  all(==(0), l.pad) || print(io, ", pad=", _maybetuple_string(l.pad))
  all(==(1), l.stride) || print(io, ", stride=", _maybetuple_string(l.stride))
  all(==(1), l.dilation) || print(io, ", dilation=", _maybetuple_string(l.dilation))
  if hasproperty(l, :groups)
    (l.groups == 1) || print(io, ", groups=", l.groups)
  end
  (l.bias === false) && print(io, ", bias=false")
end