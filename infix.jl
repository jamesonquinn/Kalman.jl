using Base.Meta

function _parseinfix(ex)
  (ex,false)
end

function parseinfix(ex)
  _parseinfix(ex)[1]
end

function _parseinfix(ex::Expr)
  if (isexpr(ex,:call)
      && ex.args[1] == :%
      && isexpr(ex.args[2],:call)
      && ex.args[2].args[1] == :%)
    return (:($(ex.args[2].args[3])($(parseinfix(ex.args[2].args[2])),
                                    $(parseinfix(ex.args[3])))),
            true)
  end
  parsedargs = map(_parseinfix,ex.args)
  (Expr(ex.head,map(first,parsedargs)...),parsedargs[1][2])
end

function addargs(ex::Expr,left,right)
  return Expr(ex.head,ex.args[1],left...,ex.args[2:end]...,right...)
end

macro %(args...)
  #println(length(a))
  if length(args) == 1
    return parseinfix(args[1])
  end
  parsedargs = map(_parseinfix,args)
  for (i, (parsedarg, changed)) in enumerate(parsedargs)
    if changed
      return(addargs(parsedarg,map(first,parsedargs[1:i-1]),
                            map(first,parsedargs[i+1:end])))
    end
  end
  throw("multiple arguments but I can't find the governing %operator%")
end
