this_file_is_not_a_function_definition = 0;
addpath(".")
global output_dir = 'build/part5/'

function y = my_heaviside(x, zero_value = 0.5)
  y = cast (x > 0, class (x));
  y (x == 0) = zero_value;
endfunction

function [Vs,D] = get_states(V, pstep, N)
  M = zeros(N, N);
  for c = 1:N
      for r = 1:N
          np = r-1-(N-1)/2;
          n  = c-1-(N-1)/2;
          if r == c
              M(r,c) = (pstep*n)^2-2*V*pstep;
          else
              M(r,c) = -2*V*sin(pstep*(np-n))./(np-n);
          end          
      end
  end
  [Vs,D] = eigenshuffle(M);
endfunction

function generate_V_graphs(V_steps, xmax, N)
  global output_dir
  pmax = 0.5*N/xmax
  pstep = 2*pmax/(N-1)
  for Vidx = 1:length(V_steps)
    V = V_steps(Vidx);

    [Vs,D] = get_states(V, pstep, N);

    X2 = ones(size(D)).';
    X1 = -X2;
    Dt = D.';
    plot([X1; X2], [Dt; Dt],'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    axis([-1.2 1.2 -6 6])
    ylabel('Energy')
    title({'Energy states of', 'a finite square well',['with V=', num2str(V,'%1.1f')]})

    set(gca,'xtick',[],'ytick',[-6:1:6])
    set(gcf, 'Position',  [1200, 400, 200, 400])
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,[output_dir, 'chart_', num2str(Vidx,'%03.f'), '.png'])
  end
endfunction

function square_well_potential(L, N, xmax)
  global output_dir
  xn = linspace(-xmax, xmax, N);
  V = 1;
  Vx = -V * my_heaviside(L-abs(xn+0.001));
  stairs(xn, Vx,'Color','b','LineWidth',2)
  set(gca,'FontSize',18)
  set(gcf, 'Position',  [1200, 400, 400, 200])
  axis([-2.5 2.5 -1.5*V 1.0*V ])
  xticks([-1, 0, 1])
  yticks([-1, 0])
  xticklabels({'-L', '0', '+L'})
  yticklabels({'-V', '0'})
  xlabel('x')
  ylabel('V')
  title('Finite square well potential')  
  print(gcf,[output_dir 'finite_square_well_potential.png'],'-dpng', '-S400,200');
endfunction

function wavefunctions(V, N, xmax, states, xaxismax, filename)
  global output_dir
  pmax = 0.5*N/xmax
  pstep = 2*pmax/(N-1)
  xn = linspace(-xmax*pi, xmax*pi, N);
  [Vs,D] = get_states(V, pstep, N);
  for idx = 1:length(states)
    phi = Vs(:,states(idx));
    pn = [-pmax:pstep:pmax];
    subplot(length(states),2,2*(idx-1)+1);
    if max(abs(real(phi))) > 0.1 * max(abs(imag(phi)))  
      plot(pn, real(phi),'Color','r','LineWidth',2)
      hold on
    end
    if max(abs(imag(phi))) > 0.1 * max(abs(real(phi)))  
      plot(pn, imag(phi),'Color','b','LineWidth',2)
    end
    axis([-5 5])
    xlabel('p')
    ylabel('Re/Im [\phi]')
    title({'\phi(p)',['with V=', num2str(V,'%1.1f')]})  
    arbitrary_factor = 0.05;
    ax = arbitrary_factor * fftshift(fft(ifftshift(phi)));  
    subplot(length(states),2,2*(idx-1)+2);
    if max(abs(real(ax))) > 0.1 * max(abs(imag(ax)))
      plot(xn, real(ax),'Color','r','LineWidth',2)
      hold on
    end
    if max(abs(imag(ax))) > 0.1 * max(abs(real(ax)))  
      plot(xn, imag(ax),'Color','b','LineWidth',2)
    end
    axis([-xaxismax xaxismax])
    xlabel('x')
    ylabel('Re/Im [\psi]')
    title({'\psi(x)',['with V=', num2str(V,'%1.1f')]})
  endfor
#  print(gcf,[output_dir filename], '-dpng', '-S600,600')
  print(gcf,[output_dir filename],
            '-dpng',
            ['-S', num2str(600,'%d'), ',', num2str(length(states) * 200,'%d')]);
endfunction

N = 401  # use odd!

#square_well_potential(m, L, N, 3*L)

L=1
xmax = 10*L
V = 1

[Vs,D] = get_states(0.5, pstep, N);
bound_state_energy = D(length(D))
wavefunctions(0.5, N, xmax, N:N, 3, 'bound_state_0.5.png')
wavefunctions(0.5, 101, 20, 1:1, 20, 'unbound_state_0.5.png')
square_well_potential(1, N, 2.5)
generate_V_graphs([0:0.1:1.5], xmax, N)
wavefunctions(3, N, xmax, N-2:N, 3, 'wavefunctions3.png')
