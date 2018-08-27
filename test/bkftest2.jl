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

using Distributions



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
    d = 2




    M = 1000
    s = 100
    if true #redo with d=8 instead of d=2
        d = 12
        bleedl = .7#1
        bleedm = .8
        bleedr = .1#.25
        hijitter = 1.
        lojitter = .25
        higap = 2
        temper = .8 #1/sqrt(2)
        basenoise = 1.
        lownoise = .16
        lowgap = 3
        noisebleed = 0.
        mciter = 6
        initialvar = 1/(1-temper)

        x0 = bkf.State(zeros(d),initialvar*eye(d))

        a = Tridiagonal(bleedl * ones(d-1),bleedm * ones(d),bleedr * ones(d-1))
        a = a * temper / (bleedl+bleedm+bleedr) #progression matrix
        #a = a * temper / det(a)^(1/d) #progression matrix - divergent

        jitvec = fill(lojitter,d)
        jitvec[1:lowgap:d] = hijitter
        g = Diagonal(jitvec)


        q = eye(d)

        f = bkf.LinearModel(a,g,q)


        noisevec = fill(basenoise,d)
        noisevec[1:lowgap:d] = lownoise
        h = Diagonal(noisevec)
        #var35 = bkf.meanvarlocs(zeros(d), inv(h), 3:5)[2]
        z = bkf.LinearObservationModel(h)
    end
    kf0 = bkf.BasicKalmanFilter(x0,f,z)

    pf = bkf.toParticleSet(kf0,M)
    pf1 = bkf.toParticleSet(kf0,1)
    pf1.particles[1:2,1] = [1, -1]
    fap = bkf.FinkelToe(pf,bkf.FinkelParams(bkf.SampleLog(3.,3.), bkf.MhSampled(1,18),1.))

    t = collect(0:2)

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
    push!(observations, bkf.Observation(pf1,1))#zeros(d)))
    print(observations[1].y, " observation\n")
    observations[1].y[1:2] = [1, -1]









    kf1 = bkf.predictupdate(kf,observations[1])
    #fra1 = bkf.FinkelParticles(fap,observations[1],1)
    fra0 = bkf.FinkelParticles(fap,observations[1],0)
    fra1k = bkf.FinkelParticles(fap,observations[1],s)
    print(fra1k.numMhAccepts/M/s/d,"\n")

    i,l,fp = [1,1,fra1k]

    #
    vis=2:min(6,d)
    cov(fra1k.tip.particles,2)[vis,vis]






    kf1.x.p[vis,vis]





    cov(fra0.tip.particles,2)[vis,vis]


    bkf.kl2(kf1.x.x,kf1.x.p,mean(fra1k.tip.particles,2),cov(fra1k.tip.particles,2))
    bkf.kl2(kf1.x.x,kf1.x.p,mean(fra1k.tip.particles,2),cov(fra1k.tip.particles,2))
    pf = bkf.toParticleSet(kf0,M^2)

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
fname = "fixed.csv"
open( fname,  "a") do outfile

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
                        "mhType",
                        "useForward",
                        "neighborhoodsize",
                        "mvlmean",
                        "mvlvar"
                        ]))
end
@load bkf

NEIGHBORHOOD_SIZE = 4
IDEAL_SAMPLES = 1000

histPerLocs = [15,30,9]
nIters = [40,0,80]
useForwards = [1.,.5,0.]
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

nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (600,600^2,div(600^2,5),10,1,20,1),
              (600,600,div(600,5)    ,40,1,20,1)]
#
nParticles = [ #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
              (80,  80^2,div(80^2,5), 4,4,20,1),
              (300,300^2,div(300^2,5),2,3,20,1),
              (800,800^2,div(800^2,5),1,2,20,1),
            (80,  80*2,div(80*2,5), 14,3,10,1),
            (300,300*2,div(300*2,5),12,2,10,1),
            (800,800*2,div(800*2,5),5,1,10,1)]
#
nParticles = [ #d,nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot,
                    #max sampType/mhType, max useForward
              #(60,80,  80^2   *10,div(80^2*2, 1), 4,2,20,1),
              #
              (36,400,400^2    ,div(400^2,5),1,1,10,1,2  ,1),
              (36,400,400^2    ,div(400^2,5),1,1,10,3  ,1,1),
              (36,400,400^2    ,div(400^2,5),1,3  ,10,1,1,1),
            (72  ,400,400^2    ,div(400^2,5),1,1,10,1,1,1),
               (36,400,40^2    ,div(40^2,5) ,1,1,10,1,1,3  ),
            (36,200,200^2      ,div(200^2,5),1,1,10,1,1,1),
            (36,200,40^2        ,div(40^2,5),5,3,10,3,2,3),
              ]

# nParticles = [(5,25,5,40,5,10,1), #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
#             (100, 100,20,10,5,10,3),
#             (500,250,10,7,3,5,2),
#             (1000,10,10,4,3,3,2)]

#Small numbers for quicker test
# nParticles = [(5,25,5,5,5,4), #nfp,npf,nfapf,reps,max nIters slot, steps
#             (100, 10000,2000,5,5,4)]
sampTypes = [bkf.SampleLog(5.,5.), bkf.SampleUniform()]
mhTypes = [bkf.MhCompromise, bkf.MhSampled]
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
mhType = mhTypes[1]
sampType = sampTypes[1]








bleedl = .7#1
bleedm = .8
bleedr = .1#.25
hijitter = 1.
lojitter = .25
higap = 2
temper = .8 #1/sqrt(2)
basenoise = 1.
lownoise = .16
lowgap = 3
noisebleed = 0.
mciter = 6
initialvar = 1/(1-temper)





for np in 1:lnParts

    d,nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc,lnParams,lnForward = nParticles[np]


    x0 = bkf.State(zeros(d),initialvar*eye(d))

    a = Tridiagonal(bleedl * ones(d-1),bleedm * ones(d),bleedr * ones(d-1))
    a = a * temper / (bleedl+bleedm+bleedr) #progression matrix
    #a = a * temper / det(a)^(1/d) #progression matrix - divergent

    jitvec = fill(lojitter,d)
    jitvec[1:lowgap:d] = hijitter
    g = Diagonal(jitvec)


    q = eye(d)

    f = bkf.LinearModel(a,g,q)


    noisevec = fill(basenoise,d)
    noisevec[1:lowgap:d] = lownoise
    r = Diagonal(noisevec)
    #var35 = bkf.meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = bkf.LinearObservationModel(Array(r))
    var35 = bkf.meanvarlocs(zeros(d), inv(r), 3:5)[2]

    kf0 = bkf.BasicKalmanFilter(x0,f,z)

    r = 1


    fpf = bkf.toParticleSet(kf0,nfp)
    pf = bkf.toParticleSet(kf0,npf)
    fapf = bkf.toFrankenSet(kf0,nfapf,NEIGHBORHOOD_SIZE)

    pf1 = bkf.toParticleSet(kf0,1)

    t = collect(0:T)


    global truth = Vector{bkf.ParticleSet}(0)#length(t))
    global observations = Vector{bkf.Observation}(0)#length(t))
    global kfs = Vector{bkf.BasicKalmanFilter}(0)#length(t))
    #global pfs = Vector{bkf.ParticleStep}(0)#length(t))
    #global faps = Vector{bkf.FrankenStep}(0)

    global kf = kf0
    push!(kfs, kf)
    global ps = bkf.ParticleStep(pf)
    #push!(pfs, ps)
    global fap = bkf.FrankenStep(fapf)
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







    r = 1

    for r in 1:reps



        open( fname,  "a") do outfile
            #
            writecsv( outfile, trsp([["truth",
                                d,#dimension
                                r,#rep number
                                1,#steps
                                "",#finkel particles
                                "",#franken particles
                                "",#part particles
                                "", #mcmc steps in finkel
                                0] #note no comma - concatenate vector
                                ["","","",""]
                                [0,#hpl
                                "none",#samptype
                                "",
                                "",
                                "",
                                "none",
                                "",
                                "",
                                bkf.meanvarlocs(truth[end],3:5)[1],0.]]))
            ##
            writecsv( outfile, trsp([["observed",
                                d,#dimension
                                r,#rep number
                                1,#steps
                                "",#finkel particles
                                "",#franken particles
                                "",#part particles
                                "", #mcmc steps in finkel
                                0]
                                 #note no comma - concatenate vector
                                ["","","",""]
                                [0,#hpl
                                "none",#samptype
                                "",
                                "",
                                "",
                                "none",
                                "",
                                "",
                                bkf.meanvarlocs(observations[end],3:5)[1],0.]]))
            ##
            writecsv( outfile, trsp([["ideal",
                                d,#dimension
                                r,#rep number
                                1,#steps
                                nfp,#finkel particles
                                nfapf,#franken particles
                                npf,#part particles
                                nIter, #mcmc steps in finkel
                                0]
                                 #note no comma - concatenate vector
                                ["","","",""]
                                [0,#hpl
                                "none",#samptype
                                "",
                                "",
                                sqe,
                                "none",
                                "",
                                "",
                                mvl[1], mvl[2]]]))
            ##
            print("\npart:")
            kls = bkf.kl2(finalDist,ps.p.particles)
            sqe = bkf.sqerr(truth[end].particles[:,1],ps)
            mvl = bkf.meanvarlocs(ps,3:5)
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
                                "none",
                                "",
                                "",
                                mvl[1], mvl[2]]]))
            print("\nfranken:")
            kls = bkf.kl2(finalDist,fap.p.particles)
            sqe = bkf.sqerr(truth[end].particles[:,1],fap)
            mvl = bkf.meanvarlocs(fap,3:5)
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
                                "none",
                                "",
                                NEIGHBORHOOD_SIZE,
                                mvl[1], mvl[2]]]))
            print("\nb14")
        end
        print("\nb15")
    ######

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
                franktime = (@timed fap = bkf.FrankenStep(fap, observations[i]))[2]
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
                writecsv( outfile, trsp([["truth",
                                    d,#dimension
                                    r,#rep number
                                    i,#steps
                                    "",#finkel particles
                                    "",#franken particles
                                    "",#part particles
                                    "", #mcmc steps in finkel
                                    0] #note no comma - concatenate vector
                                    ["","","",""]
                                    [0,#hpl
                                    "none",#samptype
                                    "",
                                    "",
                                    "",
                                    "none",
                                    "",
                                    "",
                                    bkf.meanvarlocs(truth[end],3:5)[1],0.]]))
                ##
                writecsv( outfile, trsp([["observed",
                                    d,#dimension
                                    r,#rep number
                                    i,#steps
                                    "",#finkel particles
                                    "",#franken particles
                                    "",#part particles
                                    "", #mcmc steps in finkel
                                    0]
                                     #note no comma - concatenate vector
                                    bkf.kl2(bkf.musig(finalDist)...,
                                            observations[end].y,
                                                bkf.musig(finalDist)[2])
                                    [0,#hpl
                                    "none",#samptype
                                    "",
                                    "",
                                    "",
                                    "none",
                                    "",
                                    "",
                                    bkf.meanvarlocs(observations[end],3:5)[1],
                                    var35
                                    ]]))
                ##
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
                                    "none",
                                    "",
                                    "",
                                    mvl[1], mvl[2]]]))
                #
                kls = bkf.kl2(finalDist,ps.p.particles)
                sqe = bkf.sqerr(truth[end].particles[:,1],ps)
                mvl = bkf.meanvarlocs(ps,3:5)
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
                                    "none",
                                    "",
                                    "",
                                    mvl[1], mvl[2]]]))
                #
                kls = bkf.kl2(finalDist,fap.p.particles)
                sqe = bkf.sqerr(truth[end].particles[:,1],fap)
                mvl = bkf.meanvarlocs(fap,3:5)
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
                                    "none",
                                    "",
                                    NEIGHBORHOOD_SIZE,
                                    mvl[1], mvl[2]]]))
                #
                print("\na")
            end
            print("\nb")
            print("Finkel ",r, " ... ",nfp)
            ni, nhist, nIter, i = 1,1,1,2
        end
        useforward, mhType, sampType, ni = [useForwards[1], mhTypes[1], sampTypes[1], 1]
        1
        for useForward in useForwards[1:lnForward]
            for mhType in mhTypes[1:lnParams]
                for sampType in sampTypes[1:lnParams]
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
                                                                ]))
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
end

[idealkl,finkelkl,frankenkl,partkl

]




#....
