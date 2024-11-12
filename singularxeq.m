%%
% saddle node of xeq (snxeq) should hold two cond.
% f(snxeq)=0 % right hand side of ODE
% f'(snxeq)=0
%%
function [iv,res]=singularxeq(~,iv,u)
 
z=functionofu(u')';
res=[
    z(1)-(u(iv.G)/u(iv.jnew))*u(iv.y1)-u(iv.eqb) % right hande side and left hande side for y1xeq formula  
    z(2)- (u(iv.G)/u(iv.jnew))-u(iv.eqbd)];             % diff of y1xeq formula  
end 


