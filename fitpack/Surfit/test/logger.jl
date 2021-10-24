using Logging
using LoggingExtras
function getLogger(;include=[], modulename=nothing)
    filter = include
    function smooth_filter(log_args)
        # only log with tag=smooth
         for (k,v) in log_args.kwargs
             if k==:tag && v == "smooth"
                 return true
             end
         end
         return false
     end



    function prefixFilter(log_args)
        # always return warn and error
        if log_args.level > Logging.Info
            return true
        end
    """only filter tags"""
       if isnothing(modulename) || endswith(string(log_args._module) ,modulename)
           for (k,v) in log_args.kwargs
               # func = "generateSwap
               #println("k:'$(k)' [$(eltype(k))],  v: '$(v)' [$(eltype(v))]")
               if k==:tag && v in filter
                   return true
               end
            end
            #println("$(log_args._module)::$(log_args.kwargs)")
            return false
        else
            return false
        end
    end
    smoothLogger = ActiveFilteredLogger(smooth_filter, MinLevelLogger(FileLogger("smooth.log"), Logging.Debug))
    topLogger = ActiveFilteredLogger(prefixFilter,  MinLevelLogger(FileLogger("test.log"), Logging.Debug))
    
      return  TeeLogger(smoothLogger, topLogger)
end
