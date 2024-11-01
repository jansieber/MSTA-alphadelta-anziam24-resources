function prob=powbasic_coco(prob,iv,u)
prob=coco_add_func(prob,'pw',@pw_orbit,iv,'zero','u0',u);% add equations from @pw_orbit to prob and set initial guess 
uidx_pw=coco_get_func_data(prob,'pw','uidx');% give all indices of varaibles(param)(optional 1:16)
prob=coco_add_pars(prob,'vars',uidx_pw,fieldnames(iv).');% convert all var. to param . and assign them with  
end
