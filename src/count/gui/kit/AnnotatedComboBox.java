package count.gui.kit;
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
import java.awt.Font;
import javax.swing.JComboBox;
import javax.swing.JLabel;

/**
 * A combo box with a label
 * that changes as the selection changes in the combo.
 *
 * @author csuros
 */
public class AnnotatedComboBox<E> extends JComboBox<E>
{
    /**
     *
     *
     * @param choices list of choices for the combo box
     * @param desc associated descriptions (one for each choice); will be used for the associated label.
     */
    public AnnotatedComboBox(E[] choices, String[] desc)
    {
        super(choices);
        this.desc = desc;

        initComponents();
        initListener();
    }
    private String[] desc;

    private JLabel selected_infoL;
    private void initComponents()
    {
        this.setEditable(false);

        selected_infoL = new JLabel(desc[0]);
        Font infoF = selected_infoL.getFont();
        selected_infoL.setFont(infoF.deriveFont(Font.ITALIC).deriveFont(0.8f*infoF.getSize()));
    }

    /**
     * Sets the description of the choices.
     * 
     * @param desc same length as the combobox choices
     */
    public void setChoiceInfo(String[] desc)
    {
        this.desc = desc;
        int idx = getSelectedIndex();
        if (idx<0) // shouldn't happen
        {
            selected_infoL.setText("");
        } else
        {
            selected_infoL.setText(desc[idx]);
        }
    }

    private void initListener()
    {
        this.addActionListener(action->
            {
                int idx = getSelectedIndex();
                if (idx<0) // shouldn't happen
                {
                    selected_infoL.setText("");
                } else
                {
                    selected_infoL.setText(desc[idx]);
                }
            });
    }


    public JLabel getLabel()
    {
        return selected_infoL;
    }
}
