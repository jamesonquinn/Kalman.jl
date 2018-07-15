function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

ppath("/Users/chema/mydev/Kalman.jl/src")
ppath("/Users/chema/mydev/")
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

testbkf = false
if testbkf
    @load Gadfly



    @load bkf

    x0 = bkf.State([2.,1],[1. -.99; -.99 1])
    a = eye(2)#[1. .5; -.5 1]
    a = a / det(a)
    g = [0.05 0; 0 .05]
    q = eye(2)
    f = bkf.LinearModel(a,g,q)
    h = eye(2)
    g2 = [1. 0; 0 .8]
    z = bkf.LinearObservationModel(g2)
    kf0 = bkf.BasicKalmanFilter(x0,f,z)

    pf = bkf.toParticleSet(kf0,300)
    pf1 = bkf.toParticleSet(kf0,1)
    fap = bkf.FinkelToe(pf,bkf.FinkelParams(bkf.SampleLog(3.,3.), bkf.MhCompromise(3,5,3)))

    dt = .1
    t = collect(0:dt:1)

    truth = Vector{bkf.ParticleSet}(0)#length(t))
    observations = Vector{bkf.Observation}(0)#length(t))
    kfs = Vector{bkf.BasicKalmanFilter}(0)#length(t))
    pfs = Vector{bkf.ParticleStep}(0)#length(t))
    m = [cos(tt) for tt in t]

    kf = kf0
    push!(kfs, kf)
    ps = bkf.ParticleStep(pf)
    push!(pfs, ps)
    pf

    push!(truth, pf1)
    push!(observations, bkf.Observation([0.,0]))









    kf1 = bkf.predictupdate(kf,observations[1])
    fra1 = bkf.FinkelParticles(fap,observations[1],1)
    fra100 = bkf.FinkelParticles(fap,observations[1],100)
    fra1k = bkf.FinkelParticles(fap,observations[1],1000)
    finalDist = bkf.toDistribution(kf1)
    rsamps = rand(finalDist,1000)
    mean(rsamps,2)
    #@test kf1.x == kf.x

    for i in 2:length(t)-1
        push!(truth, bkf.ap(truth[i-1]))
        push!(observations, bkf.Observation(truth[i],1))
        kf2 = bkf.predict(kf)
        kf = bkf.update(kf2,observations[i])
        push!(kfs, kf)
        ps = bkf.ParticleStep(ps,observations[i])
        push!(pfs, ps)
    end



    using DataFrames

    nshow = 30
    df1 = DataFrame(id=[x for x in pfs[end].r.s], time=length(t),
              x=vec(pfs[end].p.particles[1,:]), y=vec(pfs[end].p.particles[2,:]))
    #using Gadfly
    df2 = DataFrame(id=[x for x in pfs[end].r.s], time=length(t)- 1,
              x=vec(pfs[end-1].p.particles[1,pfs[end].r.s]),
              y=vec(pfs[end-1].p.particles[2,pfs[end].r.s]))

    df = vcat(df1[1:nshow,:],df2[1:nshow,:])
    dfh = hcat(df1,df2)

    pp=Gadfly.plot(df,x="x",y="y",group="id",Gadfly.Geom.line,
        Gadfly.Guide.annotation(
            Gadfly.compose(Gadfly.context(),
            Gadfly.ellipse(-.25,-0.5,.05,.10), Gadfly.fill(nothing),
            Gadfly.stroke("red"))) )
    inch = Gadfly.inch
    pp
    show(pp)
    Gadfly.draw(Gadfly.SVG("/Users/chema/mydev/myplot2.svg", 4inch, 3inch), pp)

    pp2=Gadfly.plot(dfh,x="x",y="x_1",Gadfly.Geom.point)
    Gadfly.draw(Gadfly.SVG("/Users/chema/mydev/myplot3.svg", 4inch, 3inch), pp2)

    3
end








#room for atom error display









function trsp(v)
    reshape(v,(1,:))
end
open( "filtertesty.csv",  "a") do outfile

    writecsv( outfile, trsp(["model",
                        "dimension",
                        "rep",
                        "steps",
                        "nfp",#finkel particles
                        "nfapf",#franken particles
                        "npf",#part particles
                        "nIter", #mcmc steps in finkel
                        "numMhAccepts",
                        "kl","covdiv","meandiv","entropydiv",
                        "histPerLoc",#hpl
                        "sampType",#samp
                        "runtime",#time
                        "errorType",
                        "sqerr",
                        "mhType"
                        ]))
end
@load bkf

NEIGHBORHOOD_SIZE = 5
IDEAL_SAMPLES = 1000

ds = [3]
ds = [25]
d = ds[end]
ln = size(ds,1)
histPerLocs = [3,9]
nIters = [20,70,200]
nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (100, 10000,2000,3,5,10,4),
              (200, 40000,8000,3,5,10,4),
              (500,250000,50000,2,3,5,3),
              (1000,    10,  10,2,3,3,3)]
#
nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (100, 10, 20,6,4,10,3),
              (200, 40,  8,3,4,10,3),
              (500,250, 50,2,3,10,2),
              (1000,10, 10,1,3,10,2)]
#
nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (80,  80^2,div(80^2,5), 4,3,10,2),
              (300,300^2,div(300^2,5),2,2,10,2),
              (800,800^2,div(800^2,5),1,2,10,2)]
#
nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (800,800,div(800,5),10,2,10,2)]
#

# nParticles = [(5,25,5,40,5,10,1), #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
#             (100, 100,20,10,5,10,3),
#             (500,250,10,7,3,5,2),
#             (1000,10,10,4,3,3,2)]

#Small numbers for quicker test
# nParticles = [(5,25,5,5,5,4), #nfp,npf,nfapf,reps,max nIters slot, steps
#             (100, 10000,2000,5,5,4)]
sampTypes = [bkf.SampleUniform(), bkf.SampleLog(5.,5.)]
mhTypes = [bkf.MhSampled, bkf.MhCompromise]
reps = max([np[4] for np in nParticles]...)#max of reps above
lnIters = length(nIters)
# finkelmeand = zeros(lnIters,ln,reps)
# frankenmeand = zeros(lnIters,ln,reps)
# partmeand = zeros(lnIters,ln,reps)
# idealmeand = zeros(lnIters,ln,reps)
finkelkl = zeros(lnIters,ln,reps)
frankenkl = zeros(lnIters,ln,reps)
partkl = zeros(lnIters,ln,reps)
idealkl = zeros(lnIters,ln,reps)

lnParts = length(nParticles)
np = 1
nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc = nParticles[np]
width = 1
mhType = mhTypes[1]
sampType = sampTypes[1]
for mhType in mhTypes
    for sampType in sampTypes

        for np in 1:lnParts

            nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc = nParticles[np]

            for width in 1:ln

                d = ds[width]
                bleed = .25
                jitter = .1
                jitterbleed = .1 # ends up being like twice this, because hits on left and right, blech.
                temper = .85
                basenoise = 1
                lownoise = .04
                lowgap = 3
                noisebleed = .1
                mciter = 6


                x0 = bkf.State(zeros(d),eye(d))

                a = SymTridiagonal(ones(d),bleed * ones(d-1))
                a = a * temper / det(a)^(1/d)

                b = SymTridiagonal(ones(d),-jitterbleed * ones(d-1))
                b[1,1] = b[d,d] = 1 - jitterbleed^2 / (1-jitterbleed^2)
                g = inv(b)
                g = g / det(g)^(1/d)

                q = eye(d)

                f = bkf.LinearModel(a,g,q)


                h = eye(d) * basenoise
                for i in 1:lowgap:d
                    h[i,i] = lownoise
                end

                b = SymTridiagonal(ones(d),-noisebleed * ones(d-1))
                b[1,1] = b[d,d] = 1 - noisebleed^2 / (1-noisebleed^2)
                r = inv(b)
                r = r / det(r)^(1/d)

                z = bkf.LinearObservationModel(h)
                kf0 = bkf.BasicKalmanFilter(x0,f,z)

                r = 1
                for r in 1:reps


                    fpf = bkf.toParticleSet(kf0,nfp)
                    pf = bkf.toParticleSet(kf0,npf)
                    fapf = bkf.toFrankenSet(kf0,nfapf,NEIGHBORHOOD_SIZE)

                    pf1 = bkf.toParticleSet(kf0,1)

                    t = collect(0:T)

                    global truth = Vector{bkf.ParticleSet}(0)#length(t))
                    global observations = Vector{bkf.Observation}(0)#length(t))
                    global kfs = Vector{bkf.BasicKalmanFilter}(0)#length(t))
                    global pfs = Vector{bkf.ParticleStep}(0)#length(t))
                    global faps = Vector{bkf.FrankenStep}(0)

                    global kf = kf0
                    push!(kfs, kf)
                    global ps = bkf.ParticleStep(pf)
                    push!(pfs, ps)
                    global fap = bkf.FrankenStep(fapf)
                    push!(faps, fap)

                    push!(truth, pf1)
                    push!(observations, bkf.Observation(pf1,1))

                    global nIter = 0
        #####
                    global finalDist = bkf.toDistribution(kfs[end])
                    # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
                    # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
                    # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
                    # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

                    print("\nideal:")
                    global rsamps = rand(finalDist,IDEAL_SAMPLES)
                    global kls = bkf.kl2(finalDist,rsamps)
                    global sqe = bkf.sqerr(truth[end].particles[:,1],rsamps)
                    global idealkl[np,width,r] = kls[1]

                    open( "filtertesty.csv",  "a") do outfile
                        writecsv( outfile, trsp([["ideal",
                                            d,#dimension
                                            r,#rep number
                                            1,#steps
                                            nfp,#finkel particles
                                            nfapf,#franken particles
                                            npf,#part particles
                                            nIter, #mcmc steps in finkel
                                            0
                                            ] #note no comma - concatenate vector
                                            kls
                                            [0,#hpl
                                            "none",#samptype
                                            "",
                                            "",
                                            sqe,
                                            "none"]]))

                        print("\npart:")
                        kls = bkf.kl2(finalDist,pfs[end].p.particles)
                        sqe = bkf.sqerr(truth[end].particles[:,1],pfs[end].p.particles)
                        partkl[np,width,r]  = kls[1]
                        writecsv( outfile, trsp([["particle",
                                            d,#dimension
                                            r,#rep number
                                            1,#steps
                                            nfp,#finkel particles
                                            nfapf,#franken particles
                                            npf,#part particles
                                            nIter,#mcmc steps in finkel
                                            0
                                            ]
                                             #note no comma - concatenate vector
                                            kls #again no comma
                                            [0,#hpl
                                            "none",#samptype
                                            "",
                                            "",
                                            sqe,
                                            "none"]]))
                        print("\nfranken:")
                        kls = bkf.kl2(finalDist,faps[end].p.particles)
                        sqe = bkf.sqerr(truth[end].particles[:,1],faps[end].p.particles)
                        frankenkl[np,width,r]  = kls[1]
                        writecsv( outfile, trsp([["franken",
                                            d,#dimension
                                            r,#rep number
                                            1,#steps
                                            nfp,#finkel particles
                                            nfapf,#franken particles
                                            npf,#part particles
                                            nIter,#mcmc steps in finkel
                                            0
                                            ]
                                             #note no comma - concatenate vector
                                            kls
                                            [0,#hpl
                                            "none",#samptype
                                            "",
                                            "",
                                            sqe,
                                            "none"]]))
            ######

                        i = 2

                        for i in 2:length(t)-1
                            push!(truth, bkf.ap(truth[i-1]))
                            push!(observations, bkf.Observation(truth[i],1))
                            kf2 = bkf.predict(kf)
                            kf = bkf.update(kf2,observations[i])
                            push!(kfs, kf)

                            parttime = (@timed ps = bkf.ParticleStep(ps,observations[i]))[2]
                            push!(pfs, ps)
                            franktime = (@timed fap = bkf.FrankenStep(fap, observations[i]))[2]
                            push!(faps, fap)

                            dif = kfs[end].x.p - kfs[end].x.p'
                            mean(dif)
                            mean(kfs[end].x.p)
                            global finalDist = bkf.toDistribution(kfs[end])
                            # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
                            # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
                            # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
                            # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

                            print("\nideal:")
                            rsamps = rand(finalDist,IDEAL_SAMPLES)
                            kls = bkf.kl2(finalDist,rsamps)
                            sqe = bkf.sqerr(truth[i].particles[:,1],rsamps)
                            idealkl[np,width,r] = kls[1]
                            writecsv( outfile, trsp([["ideal",
                                                d,#dimension
                                                r,#rep number
                                                i,#steps
                                                nfp,#finkel particles
                                                nfapf,#franken particles
                                                npf,#part particles
                                                nIter, #mcmc steps in finkel
                                                0
                                                ] #note no comma - concatenate vector
                                                kls#again no comma
                                                [0,#hpl
                                                "none",#samp
                                                "",#time
                                                "",#error
                                                sqe,
                                                "none"]]))

                            print("\npart:")
                            kls = bkf.kl2(finalDist,pfs[end].p.particles)
                            sqe = bkf.sqerr(truth[end].particles[:,1],pfs[end].p.particles)
                            partkl[np,width,r]  = kls[1]
                            writecsv( outfile, trsp([["particle",
                                                d,#dimension
                                                r,#rep number
                                                i,#steps
                                                nfp,#finkel particles
                                                nfapf,#franken particles
                                                npf,#part particles
                                                nIter,#mcmc steps in finkel
                                                0
                                                ]
                                                 #note no comma - concatenate vector
                                                kls#again no comma
                                                [0,#hpl
                                                "none",#samp
                                                parttime,#time
                                                "",#error
                                                sqe,
                                                "none"]]))
                            print("\nfranken:")
                            kls = bkf.kl2(finalDist,faps[end].p.particles)
                            sqe = bkf.sqerr(truth[end].particles[:,1],faps[end].p.particles)
                            frankenkl[np,width,r]  = kls[1]
                            writecsv( outfile, trsp([["franken",
                                                d,#dimension
                                                r,#rep number
                                                i,#steps
                                                nfp,#finkel particles
                                                nfapf,#franken particles
                                                npf,#part particles
                                                nIter,#mcmc steps in finkel
                                                0
                                                ]
                                                 #note no comma - concatenate vector
                                                kls#again no comma
                                                [0,#hpl
                                                "none",#samp
                                                franktime,#time
                                                "",#error
                                                sqe,
                                                "none"]]))
                        end
                        print("Finkel ",r, " ... ",nfp)
                        ni, nhist, nIter, i = 1,1,1,2
                        for ni = 1:lnIters
                            nIter = nIters[ni]
                            print("nIter ",nIter)
                            for nhist = 1:lnHistPerLoc
                                histPerLoc = histPerLocs[nhist]


                                fp = bkf.FinkelToe(fpf,bkf.FinkelParams(sampType,mhType(histPerLoc))) #fparams(histPerLoc,2))
                                #sampType = "sampled..uniform"
                                global fps = Vector{bkf.AbstractFinkel}(0)#length(t))
                                push!(fps, fp)
                                print("\nFinkel:",mhType," ",sampType," ",np," ",width," ",r," ",nIter," ",i)
                                kls = bkf.kl2(finalDist,fps[1].tip.particles)
                                sqe = bkf.sqerr(truth[i].particles[:,1],fps[1].tip.particles)
                                finkelkl[np,width,r]  = kls[1]
                                for i in 2:length(t)-1

                                    global finalDist = bkf.toDistribution(kfs[i])
                                    finktime = (@timed fp = bkf.FinkelParticles(fp, observations[i], nIter))[2]
                                    push!(fps, fp)
                                    print("\nfinkel:",mhType," ",sampType," ",np," ",width," ",r," ",nIter," ",i)
                                    try
                                        kls = bkf.kl2(finalDist,fps[i].tip.particles)
                                        sqe = bkf.sqerr(truth[i].particles[:,1],fps[i].tip.particles)
                                        finkelkl[np,width,r]  = kls[1]
                                        writecsv( outfile, trsp([["finkel",
                                                            d,#dimension
                                                            r,#rep number
                                                            i,#steps
                                                            nfp,#finkel particles
                                                            nfapf,#franken particles
                                                            npf,#part particles
                                                            nIter, #mcmc steps in finkel
                                                            fps[end].numMhAccepts
                                                            ] #note no comma - concatenate vector
                                                            kls #again no comma
                                                            [histPerLoc,
                                                            string(typeof(sampType)),
                                                            finktime,
                                                            "", #error
                                                            sqe,
                                                            string(mhType)
                                                            ]
                                                            ]))
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
                                                                fps[end].numMhAccepts
                                                                ] #note no comma - concatenate vector
                                                                ["","","",""] #again no comma
                                                                [histPerLoc,
                                                                string(typeof(sampType)),
                                                                finktime,
                                                                "SingularException",
                                                                "",
                                                                string(mhType)]
                                                                ]))
                                        elseif isa(y, LinAlg.LAPACKException)
                                            writecsv( outfile, trsp([["finkel",
                                                                d,#dimension
                                                                r,#rep number
                                                                i,#steps
                                                                nfp,#finkel particles
                                                                nfapf,#franken particles
                                                                npf,#part particles
                                                                nIter, #mcmc steps in finkel
                                                                fps[end].numMhAccepts
                                                                ] #note no comma - concatenate vector
                                                                ["","","",""] #again no comma
                                                                [histPerLoc,
                                                                sampType,
                                                                finktime,
                                                                "LAPACKException",
                                                                "",
                                                                string(mhType)]
                                                                ]))
                                        elseif isa(y, DomainError)
                                            writecsv( outfile, trsp([["finkel",
                                                                d,#dimension
                                                                r,#rep number
                                                                i,#steps
                                                                nfp,#finkel particles
                                                                nfapf,#franken particles
                                                                npf,#part particles
                                                                nIter, #mcmc steps in finkel
                                                                fps[end].numMhAccepts
                                                                ] #note no comma - concatenate vector
                                                                ["","","",""] #again no comma
                                                                [histPerLoc,
                                                                string(typeof(sampType)),
                                                                finktime,
                                                                "DomainError",
                                                                "",
                                                                string(mhType)]
                                                                ]))
                                            rethrow()
                                        else
                                            rethrow()
                                        end
                                    end
                                end
                            end
                        end
                    end
                end


            end
        end

    end #Mh
end #Samp
close(outfile)

[idealkl,finkelkl,frankenkl,partkl

]



wmean = zeros(d)
for l in 1:d
    wmean[l] = fps[end].ws[l]' * fps[end].base[l,:] / sum(fps[end].ws[l])
end

print([cor(wmean,bkf.params(finalDist)[1]),
cor(vec(mean(fps[end].tip.particles,2)),bkf.params(finalDist)[1]),
cor(wmean,vec(mean(fps[end].tip.particles,2))),
cor(wmean,observations[end].y),
cor(observations[end].y,bkf.params(finalDist)[1])
])


#....
