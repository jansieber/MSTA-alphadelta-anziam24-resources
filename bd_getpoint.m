function u=bd_getpoint(run,type,varnames)
if ischar(run)
    run=coco_bd_read(run);
end
inirow=find(strcmp(coco_bd_col(run,'TYPE'),type));
uall=coco_bd_col(run,varnames);
u=uall(:,inirow);
end
