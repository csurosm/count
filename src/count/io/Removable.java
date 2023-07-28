package count.io;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Interface for anything that can be removed.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public interface Removable 
{
    /**
     * Prepares the removal of this object (useful e.g., if there are some associated threads that need to be stopped)
     *
     * @return whether the operation was successful
     */
    public boolean remove();
    
}
