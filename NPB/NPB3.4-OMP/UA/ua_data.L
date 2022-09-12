NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
FORMAT LIST

 LINE   LOOP      STATEMENT

   295:                 subroutine alloc_space
   296:           
   297:           c---------------------------------------------------------------------
   298:           c---------------------------------------------------------------------
   299:           
   300:           c---------------------------------------------------------------------
   301:           c allocate space dynamically for data arrays
   302:           c---------------------------------------------------------------------
   303:           
   304:                 use ua_data
   305:                 implicit none
   306:           
   307:                 integer ios
   308:           
   309:           
   310:                 allocate (
   311:                &        ta1  (lx1,lx1,lx1,lelt), ta2   (lx1,lx1,lx1,lelt),
   312:                &        trhs (lx1,lx1,lx1,lelt), t     (lx1,lx1,lx1,lelt),
   313:                &        tmult(lx1,lx1,lx1,lelt), dpcelm(lx1,lx1,lx1,lelt),
   314:                &        pdiff(lx1,lx1,lx1,lelt), pdiffp(lx1,lx1,lx1,lelt),
   315:                >        stat = ios)
   316:           
   317:                 if (ios .eq. 0) allocate (
   318:                &        umor(lmor), tmmor(lmor),
   319:                &        rmor(lmor), dpcmor (lmor), pmorx(lmor), ppmor(lmor),
   320:                &        mormult(lmor), tmort(lmor),
   321:                >        stat = ios)
   322:           
   323:                 if (ios .eq. 0) allocate (
   324:                &        idmo    (lx1,lx1,lnje,lnje,nsides,lelt),
   325:                &        idel    (lx1,lx1,          nsides,lelt),
   326:                &        sje     (2,2,              nsides,lelt),
   327:                &        sje_new (2,2,              nsides,lelt),
   328:                &        ijel    (2,                nsides,lelt),
   329:                &        ijel_new(2,                nsides,lelt),
   330:                &        cbc     (                  nsides,lelt),
   331:                &        cbc_new (                  nsides,lelt),
   332:                &        vassign (8,lelt),       emo(2,8,8*lelt),
   333:                &        nemo    (8*lelt),
   334:                &        diagn   (2,12,lelt),
   335:                >        stat = ios)
   336:           
   337:                 if (ios .eq. 0) allocate (
   338:                &        tree   (lelt), mt_to_id    (lelt),
   339:                &        newc   (lelt), mt_to_id_old(lelt),
   340:                &        newi   (lelt), id_to_mt    (lelt),
   341:                &        newe   (lelt), ref_front_id(lelt),
   342:                &        front  (lelt), action      (lelt),
   343:                &        ich    (lelt), size_e      (lelt),
   344:                &        treenew(lelt),
   345:                >        stat = ios)
   346:           
   347:                 if (ios .eq. 0) allocate (
   348:                &        ifpcmor  (8* lelt),
   349:                &        eassign  (12,lelt),  if_1_edge(12,lelt),
   350:                &        ncon_edge(12,lelt),
   351:                &        skip (lelt), ifcoa (lelt), ifcoa_id(lelt),
   352:                &        fassign(nsides,lelt), edgevis(4,nsides,lelt),
   353:                >        stat = ios)
   354:           
   355:           c.....coordinate of element vertices
   356:                 if (ios .eq. 0) allocate (
   357:                &        xc    (8,lelt),yc    (8,lelt),zc    (8,lelt),
   358:                &        xc_new(8,lelt),yc_new(8,lelt),zc_new(8,lelt),
   359:                >        stat = ios)
   360:           
   361:                 if (ios .ne. 0) then
   362:                    write(*,*) 'Error encountered in allocating space'
   363:                    stop
   364:                 endif
   365:           
   366:                 return
   367:                 end


