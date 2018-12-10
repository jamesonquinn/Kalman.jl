function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

cd(joinpath(splitdir(@__FILE__)[1],"../output"))

ppath("/Users/chema/Dropbox/Kalman.jl/src")
ppath("/Users/chema/Dropbox/")
ppath("/Users/chema/mydev/Gadfly.jl/src")

macro load(pkg)
  quote
    if isdefined(Symbol($(string(pkg))))
      #try
      #  reload(ucfirst(string(pkg)))
      #catch
        reload($(string(pkg)))
      #end
    else
      import $pkg
    end
  end
end

using Distributions

@load bkf




#room for atom error display








valfname = "lr96fixedest10.2.csv"

function trsp(v)
    reshape(v,(1,:))
end
@load bkf

NEIGHBORHOOD_SIZE = 4
IDEAL_SAMPLES = 1000

histPerLocs = [45,30,9]
nIters = [160,0,80,20,160]
useForwards = [1.,.5,0.]
#
if false #false for quickie test
    nParticles = [ #d,nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot,
                    #max sampType/mhType, max useForward
              #(60,80,  80^2   *10,div(80^2*2, 1), 4,2,20,1),
              #
              (9,800,40^2      ,div(800^2,5),20,1,20,1,1,1),
              ]
else
# nParticles = [(5,25,5,40,5,10,1), #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
#             (100, 100,20,10,5,10,3),
#             (500,250,10,7,3,5,2),
#             (1000,10,10,4,3,3,2)]
    nParticles = [ #d,nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot,
                    #max sampType/mhType, max useForward
              #(60,80,  80^2   *10,div(80^2*2, 1), 4,2,20,1),
              #
              (9,80,4^2      ,div(80^2,5),20,1,20,1,1,1),
              ]
end

#Small numbers for quicker test
# nParticles = [(5,25,5,5,5,4), #nfp,npf,nfapf,reps,max nIters slot, steps
#             (100, 10000,2000,5,5,4)]
sampTypes = [bkf.SampleLog(5.,5.), bkf.SampleUniform()]
mhTypes = [bkf.MhSampled, bkf.MhCompromise]
sampTypes = [bkf.SampleUniform(), bkf.SampleLog(5.,5.)]
mhTypes = [bkf.MhSampled]
#mhTypes = []
reps = max([np[4] for np in nParticles]...)#max of reps above
lnIters = length(nIters)
# finkelmeand = zeros(lnIters,ln,reps)
# frankenmeand = zeros(lnIters,ln,reps)
# partmeand = zeros(lnIters,ln,reps)
# idealmeand = zeros(lnIters,ln,reps)

lnParts = length(nParticles)
finkelkl = zeros(lnParts,1,reps)
frankenkl = zeros(lnParts,1,reps)
partkl = zeros(lnParts,1,reps)
idealkl = zeros(lnParts,1,reps)
np = 1
nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc = nParticles[np]
width = 1
try
    mhType = mhTypes[1]
catch
    mhTypes = bkf.MhSampled
end
sampType = sampTypes[1]







forcingF = 8.
#d is set by nParticles above
timeStep = 0.01
numSteps = 20
processNoiseVar = 0.0001 #Is this good? Needs testing.
measurementNoiseVar = 0.0001 #Again, ???
initialvar = 0.0001 #leaves room for early progress



basenoise = .05
highnoise = .5
highgap = 5



np = 1

#for np in 1:lnParts

    d,nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc,lnParams,lnForward = nParticles[np]


    x0 = bkf.State(ones(d)*forcingF,initialvar*eye(d))

    g = eye(d) * processNoiseVar

    q = eye(d)

    f = bkf.LorenzModel(forcingF,d,timeStep,numSteps,g,q)

    if false
        #Uniform noise
        r = eye(d)*measurementNoiseVar
    else
        #low noise, with exceptions
        noisevec = fill(basenoise,d)
        noisevec[1:highgap:d] = highnoise
        r = Diagonal(noisevec)
    end

    #var35 = bkf.meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = bkf.LinearObservationModel(Array(r))
    var35 = bkf.meanvarlocs(zeros(d), inv(r), 3:5)[2]

    kf0 = bkf.BasicLorenzFilter(x0,f,z)

    r = 1


    fpf = bkf.ParticleSet(kf0,nfp)
    pf = bkf.ParticleSet(kf0,npf)
    fapf = bkf.toFrankenSet(kf0,nfapf,NEIGHBORHOOD_SIZE)

    pf1 = bkf.ParticleSet(kf0,1)

    t = collect(0:T)


    global truth = Vector{bkf.ParticleSet}(0)#length(t))
    global observations = Vector{bkf.Observation}(0)#length(t))
    global kfs = Vector{bkf.BasicLorenzFilter}(0)#length(t))
    #global pfs = Vector{bkf.ParticleStep}(0)#length(t))
    #global faps = Vector{bkf.FrankenStep}(0)

    global kf = kf0
    push!(kfs, kf)
    global ps = bkf.ParticleStep(pf)
    #push!(pfs, ps)
    global fap = bkf.FrankenStep(fapf)
    print(fap.p.particles[1:2,1:2], " fap\n")
    print(bkf.musig(fap)[2][1:2,1:2], " fapμΣ\n")
    global fapSamps = bkf.FrankenStep(fap) #null resample
    #push!(faps, fap)

    push!(truth, pf1)
    push!(observations, bkf.Observation(pf1,1))

    global nIter = 0
    #####
    global finalDist = bkf.toDistribution(kfs[end])
    pfinal = params(finalDist)
    # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
    # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
    # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
    # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

    print("\nideal:")
    global rsamps = rand(finalDist,IDEAL_SAMPLES)

    global kls = bkf.kl2(finalDist,rsamps)
    global sqe = bkf.sqerr(truth[end].particles[:,1],rsamps)
    global mvl = bkf.meanvarlocs(kfs[end],3:5)
    global idealkl[np,width,r] = kls[1]

    for i in 2:1000


        push!(truth, bkf.ap(truth[i-1]))
        print(truth[i].particles,"\n")
        print("   ",bkf.propagateUncertainty(truth[i])[1:2,1:5],"\n")
    end


if false



    r = 1

    #for r in 1:reps



        #open( fname,  "a") do outfile
            #
            ##
            print("\npart:")
            kls = bkf.kl2(params(finalDist)...,bkf.musig(ps)...)
            sqe = bkf.sqerr(truth[end].particles[:,1],ps)
            mvl = bkf.meanvarlocs(ps,3:5)
            partkl[np,width,r]  = kls[1]
            print("\nfranken:")
            print(bkf.musig(fap)[2][3:5,3:5]," fmat\n")
            klsideal = bkf.kl2(params(finalDist)...,bkf.musig(fap)...)
            kls = bkf.kl2(finalDist,fapSamps.p.particles)
            sqe = bkf.sqerr(truth[end].particles[:,1],fapSamps.p.particles)
            mvl = bkf.meanvarlocs(fap,3:5)
            frankenkl[np,width,r]  = kls[1]
            print("\nb14")
        #end
        print("\nb15")

        i = 2

        for i in 2:length(t)-1

            open( fname,  "a") do outfile
                push!(truth, bkf.ap(truth[i-1]))
                push!(observations, bkf.Observation(truth[i],1))
                print("\nideal:")
                kf2 = bkf.predict(kf)
                kf = bkf.update(kf2,observations[i])
                push!(kfs, kf)

                print("\npart:")
                parttime = (@timed ps = bkf.ParticleStep(ps,observations[i]))[2]
                #push!(pfs, ps)
                print("\nfranken:")
                franktime = (@timed begin
                    fap = bkf.FrankenStep(fapSamps, observations[i])
                    fapSamps = bkf.FrankenStep(fap)
                end)[2]
                #push!(faps, fap)

                print("\nwriting:")
                dif = kfs[end].x.p - kfs[end].x.p'
                mean(dif)
                mean(kfs[end].x.p)
                global finalDist = bkf.toDistribution(kfs[end])
                # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
                # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
                # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
                # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

                rsamps = rand(finalDist,IDEAL_SAMPLES)
                kls = bkf.kl2(finalDist,rsamps)
                sqe = bkf.sqerr(truth[i].particles[:,1],rsamps)
                mvl = bkf.meanvarlocs(kfs[end],3:5)
                idealkl[np,width,r] = kls[1]
                print("meanvarlocs(truth  ",bkf.meanvarlocs(truth[end],3:5)[1],"\n")
                ##
                print("meanvarlocs(observations[end]  ",bkf.meanvarlocs(observations[end],3:5)[1],"\n")
                kls = bkf.kl2(params(finalDist)...,bkf.musig(ps)...)
                sqe = bkf.sqerr(truth[end].particles[:,1],ps)
                mvl = bkf.meanvarlocs(ps,3:5)
                partkl[np,width,r]  = kls[1]


                kls = bkf.kl2(params(finalDist)...,bkf.musig(fap)...)
                sqe = bkf.sqerr(truth[end].particles[:,1],fap)
                mvl = bkf.meanvarlocs(fap,3:5)
                frankenkl[np,width,r]  = kls[1]
                print("\na")
            end
            print("\nb")
            print("Finkel ",r, " ... ",nfp)
            ni, nhist, nIter, i = 1,1,1,2
        end
        useforward, mhType, sampType, ni = [useForwards[1], mhTypes[1], sampTypes[1], 1]
        1
        if false #test finkel
            for useForward in useForwards[1:lnForward]
                for mhType in mhTypes[1:min(length(mhTypes),lnParams)]
                    for sampType in sampTypes[1:min(length(sampTypes),lnParams)]
                        for ni = 1:lnIters
                            nIter = nIters[ni]
                            print("nIter ",nIter)
                            for nhist = 1:lnHistPerLoc
                                histPerLoc = histPerLocs[nhist]


                                global fp = bkf.FinkelToe(fpf,bkf.FinkelParams(sampType,mhType(histPerLoc),useForward)) #fparams(histPerLoc,2))
                                #sampType = "sampled..uniform"
                                #fps = Vector{bkf.AbstractFinkel}(0)#length(t))
                                #push!(fps, fp)
                                print("\nFinkel:",mhType," ",sampType," ",np," ",width," ",r," ",nIter," ",i," ",histPerLoc," \n")
                                kls = bkf.kl2(finalDist,fp.tip.particles)
                                sqe = bkf.sqerr(truth[i].particles[:,1],fp.tip.particles)
                                finkelkl[np,width,r]  = kls[1]
                                for i in 2:length(t)-1

                                    open( fname,  "a") do outfile
                                        finalDist = bkf.toDistribution(kfs[i])
                                        finktime = (@timed fp = bkf.FinkelParticles(fp, observations[i], nIter))[2]
                                        #push!(fps, fp)
                                        print("\nfinkel:",mhType," ",sampType," ",np," ",width," ",r," ",nIter," ",i," ",histPerLoc," \n")
                                        try
                                            kls = bkf.kl2(finalDist,fp.tip.particles)
                                            sqe = bkf.sqerr(truth[i].particles[:,1],fp.tip.particles)
                                            mvl = bkf.meanvarlocs(fp,3:5)
                                            finkelkl[np,width,r]  = kls[1]
                                            writecsv( outfile, trsp([["finkel",
                                                                d,#dimension
                                                                r,#rep number
                                                                i,#steps
                                                                nfp,#finkel particles
                                                                nfapf,#franken particles
                                                                npf,#part particles
                                                                nIter, #mcmc steps in finkel
                                                                fp.numMhAccepts
                                                                ] #note no comma - concatenate vector
                                                                kls #again no comma
                                                                [histPerLoc,
                                                                string(typeof(sampType)),
                                                                finktime,
                                                                "", #error
                                                                sqe,
                                                                string(mhType),
                                                                string(useForward),
                                                                "",
                                                                mvl[1], mvl[2]
                                                                ]
                                                                ["","","",""]]))
                                        catch y
                                            if isa(y, LinAlg.SingularException)
                                                writecsv( outfile, trsp([["finkel",
                                                                    d,#dimension
                                                                    r,#rep number
                                                                    i,#steps
                                                                    nfp,#finkel particles
                                                                    nfapf,#franken particles
                                                                    npf,#part particles
                                                                    nIter, #mcmc steps in finkel
                                                                    fp.numMhAccepts
                                                                    ] #note no comma - concatenate vector
                                                                    ["","","",""] #again no comma
                                                                    [histPerLoc,
                                                                    string(typeof(sampType)),
                                                                    finktime,
                                                                    "SingularException",
                                                                    "",
                                                                    string(mhType),
                                                                    string(useForward),
                                                                    ""]

                                                                    ["","","",""]]))
                                            elseif isa(y, LinAlg.LAPACKException)
                                                writecsv( outfile, trsp([["finkel",
                                                                    d,#dimension
                                                                    r,#rep number
                                                                    i,#steps
                                                                    nfp,#finkel particles
                                                                    nfapf,#franken particles
                                                                    npf,#part particles
                                                                    nIter, #mcmc steps in finkel
                                                                    fp.numMhAccepts
                                                                    ] #note no comma - concatenate vector
                                                                    ["","","",""] #again no comma
                                                                    [histPerLoc,
                                                                    sampType,
                                                                    finktime,
                                                                    "LAPACKException",
                                                                    "",
                                                                    string(mhType),
                                                                    string(useForward),
                                                                    ""]

                                                                    ["","","",""]]))
                                            elseif isa(y, DomainError)
                                                writecsv( outfile, trsp([["finkel",
                                                                    d,#dimension
                                                                    r,#rep number
                                                                    i,#steps
                                                                    nfp,#finkel particles
                                                                    nfapf,#franken particles
                                                                    npf,#part particles
                                                                    nIter, #mcmc steps in finkel
                                                                    fp.numMhAccepts
                                                                    ] #note no comma - concatenate vector
                                                                    ["","","",""] #again no comma
                                                                    [histPerLoc,
                                                                    string(typeof(sampType)),
                                                                    finktime,
                                                                    "DomainError",
                                                                    "",
                                                                    string(mhType),
                                                                    string(useForward),
                                                                    ""]

                                                                    ["","","",""]]))
                                                rethrow()
                                            else
                                                rethrow()
                                            end
                                        end
                                    end
                                end
                            end

                        end #Mh
                    end #Samp
                end
            end
        end


    #end
end







#room for atom error display
