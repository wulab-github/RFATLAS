      program mainRandForest
      implicit none

      integer ID_num
      parameter (ID_num=572)

      integer restype
      parameter (restype=20)

      integer CGrestype
      parameter (CGrestype=7)
      integer InteractionType
      parameter (InteractionType=restype*(restype+1)/2)
      integer CGInteractionType
      parameter (CGInteractionType=CGrestype*(CGrestype+1)/2)
	  
      integer mdim,ntrain,nclass,ntest

      parameter(mdim=CGrestype*(CGrestype+1)*2)
      parameter(ntrain=ID_num-1,nclass=2,ntest=1)


      real*8 Interface_SeqInd(ID_num,mdim)
      real*8 InterfaceNorm

      real*8 x(mdim,ntrain),xts(mdim,ntest)
      integer cl(ntrain),clts(ntest)
      integer jests(ntest)
      real*8 qts(nclass,ntest)

      integer labelts,labeltr,nrnn
      integer nscale,nprot,imp,interact
      integer nsample,nrnodes,mimp,near,nprox
      integer ifprot,ifscale,iftest,mdim0,ntest0,nprot0,nscale0

      parameter(labelts=1,labeltr=1,imp=1,interact=0,
     &     nprox=1,nrnn=ntrain,nscale=0,nprot=0)
      parameter(
     &     nsample=(2-labeltr)*ntrain,
     &     nrnodes=2*nsample+1,
     &     mimp=imp*(mdim-1)+1,
     &     ifprot=nprot/(nprot-.1),
     &     ifscale=nscale/(nscale-.1),
     &     iftest=ntest/(ntest-.1),
     &     nprot0=(1-ifprot)+nprot,
     &     nscale0=(1-ifscale)+nscale,
     &     ntest0=(1-iftest)+ntest,
     &     mdim0=interact*(mdim-1)+1,
     &     near=nprox*(nsample-1)+1)
      
      real*8 Sensitivity_RF,Specificity_RF,Precision_RF,Accuracy_RF
      integer TP_RF,TN_RF,FP_RF,FN_RF
      integer Neg_Acc,Pos_Acc,Bind_Ind

      integer id,index,index2
      integer i,j,ii,jj,k
      integer Interface_Contact_Map_core(InteractionType,ID_num)
      integer Interface_Contact_Map_MH(InteractionType,ID_num)
      integer Interface_Contact_Map_PG(InteractionType,ID_num)
      integer Interface_Contact_Map_MG(InteractionType,ID_num)
      character*4 pdbid(ID_num)
      character*2 WMtype(ID_num)
      real*8 affinity(ID_num)
      integer contact_num(ID_num)
      integer CG_typeindex(restype)
      integer AA_CG_mapping(InteractionType)	  
	  
cccccccccccccccc


      open(unit=11,file='AminoAcids_7type.dat',
     &     status='old')
      do i=1,restype
         read(11,2100) CG_typeindex(i)
      enddo
 2100 format(3x,I3)
      close(11)

      index=0
      do i=1,CGrestype
         do j=1,i
            index=index+1
            index2=0
            do ii=1,restype
               do jj=1,ii
                  index2=index2+1
                  if(CG_typeindex(ii).ge.CG_typeindex(jj))then
                     if((CG_typeindex(ii).eq.i).AND.
     &                    (CG_typeindex(jj).eq.j))then
                        AA_CG_mapping(index2)=index
                     endif
                  elseif(CG_typeindex(ii).lt.CG_typeindex(jj))then
                     if((CG_typeindex(jj).eq.i).AND.
     &                    (CG_typeindex(ii).eq.j))then
                        AA_CG_mapping(index2)=index
                     endif
                  endif
               enddo
            enddo

         enddo
      enddo

      open(unit=10,file=
     &        'ATLASAff_InterfaceStat_4Region_02252021.dat',
     &        status='old')

      do id=1,ID_num

         contact_num(id)=0

         read(10,2200) pdbid(id),affinity(id),WMtype(id)   
 2200    format(A4,1x,F6.2,1x,A2)
         
         do i=1,InteractionType
            read(10,2300)
     &           Interface_Contact_Map_core(i,id),
     &           Interface_Contact_Map_MH(i,id),
     &           Interface_Contact_Map_PG(i,id),
     &           Interface_Contact_Map_MG(i,id)
         enddo

 2300    format(8x,4I10)
         
      enddo
      close(10)
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	  
c>>   construct the input vector Interface_SeqInd(ID_num,mdim)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  
      do id=1,ID_num
         do j=1,mdim
            Interface_SeqInd(id,j)=0
         enddo
      enddo	  
	  
      do id=1,ID_num
         do j=1,InteractionType
            Interface_SeqInd(id,AA_CG_mapping(j))=
     &           Interface_SeqInd(id,AA_CG_mapping(j))
     &           +Interface_Contact_Map_core(j,id)
         enddo
         do j=1,InteractionType
            Interface_SeqInd(id,AA_CG_mapping(j)
     &           +CGInteractionType)=
     &           Interface_SeqInd(id,AA_CG_mapping(j)
     &           +CGInteractionType)
     &           +Interface_Contact_Map_MH(j,id)
         enddo
         do j=1,InteractionType
            Interface_SeqInd(id,AA_CG_mapping(j)
     &           +2*CGInteractionType)=
     &           Interface_SeqInd(id,AA_CG_mapping(j)
     &           +2*CGInteractionType)
     &           +Interface_Contact_Map_PG(j,id)
         enddo
         do j=1,InteractionType
            Interface_SeqInd(id,AA_CG_mapping(j)
     &           +3*CGInteractionType)=
     &           Interface_SeqInd(id,AA_CG_mapping(j)
     &           +3*CGInteractionType)
     &           +Interface_Contact_Map_MG(j,id)
         enddo	 		 
      enddo	  	  
	  
      do id=1,ID_num
         InterfaceNorm=0
         do j=1,mdim
            InterfaceNorm=InterfaceNorm+Interface_SeqInd(id,j)
         enddo
         do j=1,mdim
            Interface_SeqInd(id,j)=Interface_SeqInd(id,j)/InterfaceNorm
         enddo
      enddo	  	  
      
cccccccccccccccccccccccccccccccccccccccccccccccc	  

      TP_RF=0
      FP_RF=0
      FN_RF=0
      TN_RF=0

      print*,"Entry id   TP   TN   FP   FN"

      do id=1,ID_num

 
cccccccccc   Random Forest

         index=0
         do i=1,ID_num 
            if(i.ne.id)then
               index=index+1
               do j=1,mdim
                  x(j,index)=Interface_SeqInd(i,j)
               enddo
               if(affinity(i).ge.-6.45)then
                  cl(index)=1
               elseif(affinity(i).lt.-6.45)then
                  cl(index)=2
               endif
            endif
         enddo

         do j=1,mdim
            xts(j,1)=Interface_SeqInd(id,j)
         enddo
         if(affinity(id).ge.-6.45)then
            clts(1)=1
         elseif(affinity(id).lt.-6.45)then
            clts(1)=2
         endif

         do i=1,ntest
            jests(i)=0
            qts(1,i)=0
            qts(2,i)=0
         enddo

         call RandForest(mdim,ntrain,nclass,
     &        ntest,x,cl,xts,clts,
     &        jests,qts,
     &        labelts,labeltr,
     &        nrnn,interact,imp,nprot,nscale,nprox,
     &        nsample,nrnodes,mimp,near,
     &        ifprot,ifscale,iftest,mdim0,ntest0,nprot0,nscale0)
         
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if((affinity(id).ge.-6.45)
     &        .AND.(jests(1).eq.1))then
            TN_RF=TN_RF+1
         elseif((affinity(id).ge.-6.45)
     &        .AND.(jests(1).eq.2))then
            FP_RF=FP_RF+1
         elseif((affinity(id).lt.-6.45)
     &        .AND.(jests(1).eq.2))then
            TP_RF=TP_RF+1
         elseif((affinity(id).lt.-6.45)
     &        .AND.(jests(1).eq.1))then
            FN_RF=FN_RF+1
         endif

         print*,id,TP_RF,TN_RF,FP_RF,FN_RF


      enddo
ccccccccccccccccccccccccccccccc

      Sensitivity_RF=real(TP_RF)/real(TP_RF+FN_RF)
      Specificity_RF=real(TN_RF)/real(TN_RF+FP_RF)
      Precision_RF=real(TP_RF)/real(TP_RF+FP_RF)
      Accuracy_RF=real(TP_RF+TN_RF)/real(TP_RF+TN_RF+FP_RF+FN_RF)


      print*,'Sensitivity: ',Sensitivity_RF
      print*,'Specificity: ',Specificity_RF
      print*,'Precision: ',precision_RF
      print*,'Accuracy: ',Accuracy_RF

ccccccccccccccccccccccccccc

      stop
      end

