using Cairo, Fontconfig, Gadfly


img = plot(Theme(background_color = "white"), Guide.annotation()
    layer(x=newgr[:,1], y=newrmf[:,1], Geom.point, Geom.line, Theme(default_color=colorant"green")),
    layer(x=newgr[:,2], y=newrmf[:,2], Geom.point, Geom.line, Theme(default_color=colorant"blue")),  
    layer(x=newgr[:,3], y=newrmf[:,3], Geom.point, Geom.line, Theme(default_color=colorant"pink")),
    layer(x=newgr[:,4], y=newrmf[:,4], Geom.point, Geom.line, Theme(default_color=colorant"yellow")),
    layer(x=newgr[:,5], y=newrmf[:,5], Geom.point, Geom.line, Theme(default_color=colorant"red")),
    layer(x=newgr[:,6], y=newrmf[:,6], Geom.point, Geom.line, Theme(default_color=colorant"purple"))    )

draw(PNG("chl_ns.png", 15cm, 9cm), img)