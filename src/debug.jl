
canDebug = true
if canDebug

    verbose = true
    function debug(a...)
        global verbose
        if verbose
            for i=a
                show(i)
                print(" ")
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
            debug(a...)
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
