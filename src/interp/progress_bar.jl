@kwdef mutable struct ProgressBar
    i::Int = 0
    total::Int = -1
    desc::String = "Progress: "
    offset::Int = 0
    enabled::Bool = true
    finished::Bool = false
    terminal::Base.Terminals.TTYTerminal = Base.Terminals.TTYTerminal(string(), stdin, stdout, stderr)
    time_start::Float64 = time()
end

progress_bar_lock = ReentrantLock()

function update!(prog::ProgressBar, i::Int)
    lock(progress_bar_lock) do
        prog.i = i
        if !prog.enabled
            return
        end
        buffer = prog.desc * "        " * string(prog.i) * " / "
        if prog.total > 0
            buffer *= string(prog.total)
        else
            buffer *= "?"
        end
        if prog.i > 0
            time_end = time()
            elapsed = time_end - prog.time_start
            buffer *= " (" * @sprintf("%.3f", elapsed) * " s)"
            if prog.finished
                buffer *= " (done)"
            end
        end
        for _ in 1:prog.offset
            if prog.i == 0
                println()
            else
                Base.Terminals.cmove_line_down(prog.terminal)
            end
        end
        Base.Terminals.clear_line(prog.terminal)
        print(buffer)
        Base.Terminals.cmove_left(prog.terminal, length(buffer))
        if prog.finished
            if prog.offset == 0
                println()
            end
        end
        for _ in 1:prog.offset
            Base.Terminals.cmove_line_up(prog.terminal)
        end
    end
end

function next!(prog::ProgressBar)
    update!(prog, prog.i + 1)
end

function finish!(prog::ProgressBar)
    prog.finished = true
    if prog.total == -1
        prog.total = prog.i
    end
    update!(prog, prog.i)
end
