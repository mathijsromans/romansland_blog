this_file_is_not_a_function_definition = 0;
addpath("..")

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

function generate_V_graphs(V_steps, pstep, pmax, N)
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
    saveas(gcf,['chart_', num2str(Vidx,'%03.f'), '.png'])
  end
endfunction

function square_well_potential(m, L, N, xmax)
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
  print(gcf,'finite_square_well_potential.png','-dpng', '-S400,200');
endfunction

function wavefunctions(V, N, pstep, pmax, xmax)
  xn = linspace(-xmax*pi, xmax*pi, N);
  [Vs,D] = get_states(V, pstep, N);
  num_states = 3;
  for idx = 1:num_states
    phi = Vs(:,length(D)-num_states+idx);
    pn = [-pmax:pstep:pmax];
    subplot(num_states,2,2*(idx-1)+1);
    if max(abs(real(phi))) > 0.1 * max(abs(imag(phi)))  
      plot(pn, real(phi),'Color','r','LineWidth',2)
      hold on
    end
    if max(abs(imag(phi))) > 0.1 * max(abs(real(phi)))  
      plot(pn, imag(phi),'Color','b','LineWidth',2)
    end
    axis([-10 10])
    xlabel('p')
    ylabel('Amplitude')
    title({'Finite square well',['with V=', num2str(V,'%1.1f'),' n=', num2str(length(D)-3+idx)]})  
    arbitrary_factor = 0.1;
    ax = arbitrary_factor * fftshift(fft(ifftshift(phi)));  
    subplot(num_states,2,2*(idx-1)+2);
    if max(abs(real(ax))) > 0.1 * max(abs(imag(ax)))
      plot(xn, real(ax),'Color','r','LineWidth',2)
      hold on
    end
    if max(abs(imag(ax))) > 0.1 * max(abs(real(ax)))  
      plot(xn, imag(ax),'Color','b','LineWidth',2)
    end
    axis([-2 2])
    xlabel('x')
    ylabel('Re/Im')
    title({'Finite square well',['with V=', num2str(V,'%1.1f'),' n=', num2str(length(D)-3+idx)]})
  endfor
  print(gcf,'wavefunctions.png','-dpng', '-S600,600');
endfunction

N = 401  # use odd!

#square_well_potential(m, L, N, 3*L)

L=1
xmax = 5*L
pmax = 0.5*N/xmax
pstep = 2*pmax/(N-1)
V = 1

#[Vs,D] = get_states(V, pstep, N);

#generate_V_graphs([0:0.1:2], pstep, pmax, N)
wavefunctions(3*V, N, pstep, pmax, xmax)