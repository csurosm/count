# creates the jar files:
# main jar includes source an javadoc; CountXXIIIslim.jar contains only the class and image files
# launch as bash make_cls_jar.bash
jar --create --file CountXXIII.jar --manifest ../manifest.mf --main-class count.von.Count -C ../build/classes count -C ../build/classes img -C .. src/count -C .. doc 
jar --create --file CountXXIIIslim.jar --manifest ../manifest.mf --main-class count.von.Count -C ../build/classes count -C ../build/classes img 
