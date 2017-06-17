# Solve the Schrodinger equation for the Stern-Gerlach experiment numerically.
# coordinate system:
#   x is the beam axis, ignored here, so this is actually S. eqn. in 2 dimensions
#   z is the transverse direction along which the magnetic field points
# Spin indices are 0 for s(z) up, 1 for s(z) down.
# Parameters can be set either by editing the code or on the command line, like this:
#   ruby sg.rb '$n = 41; packet_width=0.15; t=0.025; niter=8000; mag0=300.0; k=300.0; suffix="x"'
# Command line overrides anything set in code.
# The suffix parameter gives a way to vary output the filenames so we can run the program on two CPUs concurrently.

$PI=3.1415926535

#----------- main

def main
  # parameters:
  #   $n = number of grid points in y and z directions; should be odd
  #   packet_width = rms width of initial gaussian wave packet, in units where the whole grid is [-1,1]x[-1,1]
  #   t = time for which to run the simulation
  #   niter = number of iterations to split t into; t/niter should be small compared to inverse of field
  #           for good numerical precision
  #   mag0 = B = homogeneous part of magnetic field
  #   k = dBz/dz=-dBy/dy, with y and z measured in the same distance units as packet_width
  # units are such that hbar=1, mass=1, magnetic moment of electron=1
  # The following seems like a reasonable set of defaults to start from:
  #    $n = 41; packet_width=0.15; t=0.025; niter=2500; mag0=20.0; k=300.0

  suffix = ''

  #$n = 41; packet_width=0.15; t=0.025; niter=2500; mag0=20.0; k=300.0
  #$n = 41; packet_width=0.15; t=0.025; niter=5000; mag0=300.0; k=300.0
  # ... two clearly defined wave packets, overlapping
  #$n = 41; packet_width=0.15; t=0.04; niter=8000; mag0=300.0; k=300.0
  # ... somewhat more separated; asymmetry suggests poor numerical quality
  #$n = 21; packet_width=0.15; t=0.07; niter=20000; mag0=300.0; k=300.0
  # ... somewhat more separated; asymmetry suggests poor numerical quality
  #$n = 41; packet_width=0.15; t=0.07; niter=40000; mag0=300.0; k=300.0
  # ... this is leading to very long computation times
  $n = 41; packet_width=0.15; t=0.04; niter=20000; mag0=300.0; k=300.0

  if ARGV.length==1 then eval(ARGV[0]) end

  prec = niter/(t*(mag0.abs+k.abs))
  print "niter/(tB)=#{prec}, should be >~300\n"
  if prec<300 then exit(-1) end

  print "A=B/kr=#{mag0/(k*packet_width)} (for r=initial rms value)\n"

  print "current time               = #{Time.now}\n"
  comp_time = niter*$n*$n*(5.7e-7)
  print "estimated computation time = #{comp_time} minutes\n"
  print "ETA                        = #{Time.at(Time.now.to_i+(60*comp_time).to_i)}\n"

  init_graph = "a#{suffix}.png"
  final_graph = "b#{suffix}.png"
  print "writing to files #{init_graph} #{final_graph}\n"
  
  psi = init_psi(true,packet_width)
  #ascii_print(psi)
  graph(psi,init_graph)
  dt = t/niter # should be small compared to inverse of field
  0.upto(niter-1) { |i|
    print "\n#{i}/#{niter} " if i%(niter/10)==0
    print "." if i%(niter/100)==0
    psi=propagate(psi,dt,mag0,k,i,niter/100)
  }
  print "\n"
  #ascii_print(psi)
  graph(psi,final_graph)
end

#------------- graphics data

require 'chunky_png' # ubuntu package ruby-chunky-png

require 'hsluv' 
  # http://www.hsluv.org
  # https://github.com/hsluv/hsluv-ruby
  # sudo gem install hsluv

$gamma = 1.5 # https://en.wikipedia.org/wiki/Gamma_correction
             # Should be pretty close to 1, since hsluv is perceptual?

$default_width = 614 # 52 mm at 300 dpi
$bg_color = ChunkyPNG::Color::rgb(222,222,222)

# Define a cyclic color map using a lookup table with 12 interpolation points.
# Each entry is [h,s,v] in the hsluv system.
$color_lut = [
  [8,76,60],  # 0
  [20,77,60],  # 1
  [30,78,60],  # 2
  [57,86,69],  # 3
  [94,79,65],  # 4
  [142,67,60],  # 5
  [192,64,60],  # 6
  [235,64,60],  # 7
  [270,64,50],  # 8
  [300,64,60],  # 9
  [327,64,60],  # 10
  [348,64,60]  # 11
]

#---------------------- graphics output

def graph(psi,file)
  k = 6 # larger pixels, easier to view
  w = $n*k
  h = $n*k
  image = ChunkyPNG::Image.new(2*w,h,$bg_color)
  max = max_abs(psi)
  max2 = max*max
  0.upto(w-1) { |jy|
    0.upto(h-1) { |jz|
      0.upto(1) { |j|
        iy = jy/k
        iz = jz/k
        # z = (psi[iy][iz][0].abs2+psi[iy][iz][1].abs2)/max2 # display probability, no phase
        z = psi[iy][iz][j]/max # display wf's spin components
        c = color([z.abs,z.arg])
        if j==1 && jy==0 then c = ChunkyPNG::Color::rgb(255,255,255) end # draw a white border separating the two boxes
        image[jy+w*j,h-jz-1] = c # spin up in left box, down in right
      }
    }
  }
  image.save(file)
end


#----------- physics

def propagate(psi,dt,mag0,k,iter,debug_interval)
  # iter is iteration count, used only for debugging
  psi2 = init_psi(false)
  i = Complex(0,1)
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      0.upto(1) { |j|
        debug = false
        # debug = iter%debug_interval==0 && (iy==$n/2 && iz==$n/2)
        # Compute the Laplacian. Don't try to do anything sensible at the edges.
        laplacian = -4.0*psi[iy][iz][j]
        if iy>0 then laplacian=laplacian-psi[iy-1][iz][j] end
        if iy<$n-1 then laplacian=laplacian-psi[iy+1][iz][j] end
        if iz>0 then laplacian=laplacian-psi[iy][iz-1][j] end
        if iz<$n-1 then laplacian=laplacian-psi[iy][iz+1][j] end
        laplacian = laplacian/(dx*dx)
        kin = -0.5*laplacian # kinetic energy part of hamiltonian
        y,z = indices_to_coords(iy,iz)
        if j==0 then spin=1.0 else spin=-1.0 end
        psi_old_self = psi[iy][iz][j]
        psi_old_other = psi[iy][iz][1-j]
        h_self = (mag0+k*z)*spin*psi_old_self # diagonal term in coupling to B
        h_other = i*k*y*spin*psi_old_other    # off-diagonal term
        ham = kin+h_self+h_other
        psi2[iy][iz][j] = psi_old_self-i*ham*dt
        if debug then
          #print "psi_old_self=#{psi_old_self}, psi_old_other=#{psi_old_other}, kin=#{kin}, h_self=#{h_self}, h_other=#{h_other}\n"
          print "|psi2|^2=#{psi2[iy][iz][j].abs2}\n"
        end
      }
    }
  }
  return normalize(psi2)
end

def normalize(psi)
  psi2 = init_psi(false)
  norm = 0.0
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      0.upto(1) { |j|
        norm = norm+psi[iy][iz][j].abs2
      }
    }
  }
  norm = Math::sqrt(norm)
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      0.upto(1) { |j|
        psi2[iy][iz][j] = psi[iy][iz][j]/norm
      }
    }
  }
  return psi2
end

def max_abs(psi)
  max = 0.0
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      0.upto(1) { |j|
        z = psi[iy][iz][j].abs2
        if z>max then max=z end
      }
    }
  }
  return Math::sqrt(max)
end

def init_psi(if_fill,packet_width=0.2)
  # if_fill=true means fill it with the initial wave packet, false means leave it zero
  # for if_fill=false, can omit packet_width, which is irrelevant
  psi = Array.new($n) {Array.new($n) {Array.new(2,Complex(0,0))}}
    # wavefunction, psi[iy][iz][j], where j is a spin index
  if !if_fill then
    return psi
  end
  # This will not be normalized, but we don't care, will handle that when we display the results.
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      y,z = indices_to_coords(iy,iz)
      packet = Math::exp(-0.5*(y*y+z*z)/(packet_width*packet_width)) # Gaussian envelope
      #if y.abs<0.1 && z.abs<0.1 then print "y=#{y} z=#{z} packet=#{packet}\n" end
      0.upto(1) { |j|
        if j==0 then z=Complex(1,0) else z=Complex(1,0) end # polarized along x
        #if j==0 then z=Complex(1,0) else z=Complex(0,1) end # polarized along y
        psi[iy][iz][j] = packet*z
      }
    }
  }
  return normalize(psi)
end

def indices_to_coords(iy,iz)
  return [index_to_coords(iy),index_to_coords(iz)]
end

def index_to_coords(i) # returns a number from -1 to 1
  return 2.0*i.to_f/($n-1).to_f-1.0
end

def dx # spacing of the grid
  return 2.0/($n-1).to_f
end

def ascii_print(psi)
  print "------------------------------------------------------\n"
  0.upto($n-1) { |iy|
    0.upto($n-1) { |iz|
      print "#{print_wf(psi[iy][iz])} "
    }
    print "\n"
  }
  print "------------------------------------------------------\n"
end

def print_wf(psi)
  return "<#{print_complex(psi[0])},#{print_complex(psi[1])}>"
end

def print_complex(z)
  return "#{print_real(z.real)}+#{print_real(z.imag)}i"
end

def print_real(x)
  return "%7.4f" % x
end

#--------------- low-level color stuff

def color(z)
  z = clean_up_complex(z)
  hue,saturation,value = complex_to_hsluv(z)
  r,g,b = Hsluv::hsluv_to_rgb(hue, saturation, value)
  return ChunkyPNG::Color::rgb((r*255).to_i,(g*255).to_i,(b*255).to_i)
end

def clean_up_complex(z)
  if z[0]>=0.0 then return z end
  z[0] = -z[0]
  z[1] = put_in_range(z[1]+$PI,2.0*$PI)
  return z
end

def complex_to_hsluv(z)
  mag,arg = z
  value_scale = mag**$gamma
  h,s,v = interpolate_color_lut(arg)

  # fix glitch where blue is too bright and saturated when v is low
  blue = 266.0
  w = 6.0
  if h>blue-w && h<blue+w then
    x = 0.75 # knob to turn the correction up and down, goes 0 to 1
    u = (h-blue)/w # varies from -1 to 1 in this band
    if mag<1.0 then c=1.0-mag else c=0.0 end
    a=Math::cos(u*$PI/2.0) # runs from 0 to 1, measures how much to correct by
    v=v*(1.0-0.1*x*a*c*c)
    s=s*(1.0-1.0*x*a*c)
  end

  return [h,s,v*value_scale]
end

def interpolate_color_lut(arg)
  alpha = 12.0*put_in_range(arg/(2.0*$PI),1.0)
  if alpha==12.0 then alpha=11.999 end
  i = alpha.to_i
  j = (i+1)%12
  x = alpha-i
  h1,s1,v1 = $color_lut[i]
  h2,s2,v2 = $color_lut[j]
  if h2<h1 then h2=h2+360.0 end
  h,s,v = [interpolate(h1,h2,x),interpolate(s1,s2,x),interpolate(v1,v2,x)]
  if h>360.0 then h=h-360.0 end
  return h,s,v
end

def interpolate(a,b,x)
  return a*(1.0-x)+b*x
end

def put_in_range(x,max)
  i = (x/max).to_i
  x = x-i*max
  if x>max then x=x-max end
  if x<0 then x=x+max end
  return x
end

main
