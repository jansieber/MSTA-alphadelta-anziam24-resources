function prob=snbasic_coco(iv,u)
prob=coco_prob();
prob=coco_add_func(prob,'smallSaddlexeq',@singularxeq,iv,'zero','u0',u);
uidex=coco_get_func_data(prob,'smallSaddlexeq','uidx');
prob=coco_add_pars(prob,'vars',uidex,fieldnames(iv).');
prob = coco_add_glue(prob,'gluey0103',iv.y01,iv.y03);
end
