
canDebug = false
if canDebug

    verbose = false
    function debug(a...)
        global verbose
        if verbose
            for i=a
                print(i," ")
            end
            print("\n")
        end
    end
    function setVerbose(v=true)
        global verbose = v
    end
    debugDone = Dict{String,Bool}()
    function debugOnce(a...)
        global debugDone
        if !get(debugDone,a[1],false)
            debugDone[a[1]] = true
            for i=a
                print(i," ")
            end
            print("\n"); print("\n","""print("\n")""")
        end
    end

    function resetOnce()
        global debugDone = Dict{String,Bool}()
    end


else

    #stubs
    function debug(a...)
    end
    function debugOnce(a...)
    end
    function resetOnce()
    end

end
