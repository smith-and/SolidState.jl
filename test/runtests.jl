using SolidState
using Test

@testset "SolidState.jl" begin
    # Write your tests here.
    SolidState

    shift = [1.0,0.0,0.0]
    shifted_asd = SolidState.LayerShiftASD(ASD2,[[0.0,0.0,0.0],shift])
    hd = TightBindingDensity(shifted_asd)
    dm = DataMap(dtype,shifted_asd,hd,indices,priors,base)
    di = DataIntegral(dm)
    di(Nevals:Nevals)
    di.data[1][:]

    asd = ASD2
    asd0 = asd()
    shifted_asd = SolidState.LayerShiftASD(deepcopy(asd0),[[0.0,0.0,0.0],[1.0,0.1,0.0]])

    shifted_asd["sk"]==asd0["sk"]
    shifted_asd["blv"]==asd0["blv"]
    shifted_asd["lbase"]==asd0["lbase"]
    shifted_asd["cutoff"]==asd0["cutoff"]
    shifted_asd["scales"]==asd0["scales"]
    shifted_asd["filling"]==asd0["filling"]
    shifted_asd["regulator"]==asd0["regulator"]
    shifted_asd["sites"]==asd0["sites"]

    shifted_asd["sites"][1]==asd0["sites"][1]
    shifted_asd["sites"][2]==asd0["sites"][2]
    shifted_asd["sites"][3]==asd0["sites"][3]
    shifted_asd["sites"][4]==asd0["sites"][4]

    (shifted_asd["sites"][3][3].-asd0["sites"][3][3])==[1.0,0.0,0.0]
    
end
