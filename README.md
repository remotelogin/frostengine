# frostengine!

![frostengine_showcase_image](https://raw.githubusercontent.com/remotelogin/frostengine/main/.github_files/frostengine_showcase.gif)

---

 Its an engine! It works, its single threaded, but still performs pretty well ngl.
 
 its written in plain cpp20 and x11 to draw indiviual pixels to a window. the rest is entirely made from scratch.

 also no, multi threading it would not rly increase performance, cause im heavily limited by x11 being incredibly slow at drawing lines (see the performance statistics in the renderer's top left)
 
 # features
 
  - obj file loading to visualize models ( big models require a fast clock speed cpu. otherwise rendering slows down )
  - double buffering
  - abstracted and easy to use model loading and positioning (oop!!!! omg i love oop!!!)
  - 3 dof model transformation (per model)
  - 3 dof camera movement and "jet" style camera
  - *VERY* simple face shading
  
# to do:

 - texture rendering
 - abstract everything
 - advanced shading / lighting engine
 
# to do (far in the future lol):
 
 - ray tracing xdddd
 - lighting engine
 - animation loader
 - physics / collission engine
