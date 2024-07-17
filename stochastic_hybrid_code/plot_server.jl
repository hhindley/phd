using WGLMakie
using DataFrames
using HTTP
using WebIO

# Function to create the plot
function create_plot(data)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Example Plot")
    lines!(ax, data.x, data.y)
    return fig
end

# Example data
df = DataFrame(x = 1:1000, y = rand(1000))

# Create the plot
fig = create_plot(df)

# Serve the interactive plot
function handle_request(req)
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Interactive Plot</title>
        <script src="https://cdn.jsdelivr.net/npm/webio@0.8.0/dist/webio.min.js"></script>
    </head>
    <body>
        <div id="plot">
            $(WebIO.render($(fig)))
        </div>
        <script>
            document.addEventListener("DOMContentLoaded", function() {
                var ws_url = "ws://" + window.location.host + "/__webio_ws__";
                var webio = new WebIO.Endpoint({url: ws_url});
                WebIO.mount(webio, document.getElementById("plot"));
            });
        </script>
    </body>
    </html>
    """
    return HTTP.Response(200, html)
end

HTTP.serve(handle_request, "0.0.0.0", 8080)
